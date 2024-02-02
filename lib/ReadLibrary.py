#!/usr/bin/env python
import sys
import subprocess
import argparse
import gzip
import bz2
import os
import os.path
import re
import glob
import shutil
import copy
from time import time, localtime, strftime, sleep

def inferPlatform(read_id, maxReadLength, avgReadQuality):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    if maxReadLength < ReadLibrary.MAX_SHORT_READ_LENGTH:
        # example illumina read id
        #@D00553R:173:HG53VBCXY:2:1101:1235:2074 1:N:0:ACAGTGAT
        if read_id.startswith(">"):
            return "fasta"
        parts = read_id.split(":")
        if re.match(r"@[A-Z]\S+:\d+:\S+:\d+:\d+:\d+:\d+ \S+:\S+:\S+:\S+$", read_id):
            return "illumina" # newer illumina
        if re.match(r"@\S+:\S+:\S+:\S+:\S+#\S+/\S+$", read_id):
            return "illumina" # older illumina
        if re.match(r"@[^:]+:[^:]+:[^:]+$", read_id):
            return "iontorrent" # 
        if re.match(r"@[SED]RR\d+\.\d+", read_id):
            return "illumina" # default short fastq type 
        if len(parts) > 4:
            return "illumina"
        if len(parts) == 3:
            return "iontorrent"
    else: # one of the long read types (pacbio or nanopore)
        if read_id.startswith(">"):
            return "nanopore"
        if avgReadQuality > 11 and avgReadQuality < 31:
            sys.stderr.write("inferring platform is nanopore from quality: {:.4f}\n".format(avgReadQuality))
            return "nanopore"
        else:
            sys.stderr.write("inferring platform is pacbio from quality: {:.4f}\n".format(avgReadQuality))
            return "pacbio"
        # todo: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
        dash_delimited_fields = read_id.split("-")
        if len(dash_delimited_fields) == 5:
            return "nanopore" # based on user data
        if re.match(r"@\S+/\S+/\S+_\S+$", read_id): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
            return "pacbio" # 
        if re.match(r"@\S+/\S+$", read_id): #@<MovieName>/<ZMW_number> 
            return "pacbio" # 
    #@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
        if re.match(r"@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", read_id): #based on one example, need to test more 
            return "nanopore" # 
        return "pacbio" # default long fastq type
    # if we get here, we failed to recognize what read type it is, default to illumina
    sys.stderr.write("inferPlatform defaulting to 'illumina'\n")
    return "illumina"

def findSingleDifference(s1, s2):
    # if two strings differ in only a single contiguous region, return the start and end of region, else return None
    if len(s1) != len(s2):
        return None
    start = None
    end = None
    i = 0
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 != c2:
            if end:
                return None
            if not start:
                start = i
        elif start and not end:
            end = i
    if start and not end:
        end = i+1
        return (start, end)

class ReadLibrary:
    # representation of read set, with specific versions included within (e.g., trimmed version)
    MEMORY = '5gb'
    MAX_SHORT_READ_LENGTH = 999
    NUM_THREADS = 4
    MAX_BASES=1e9
    LOG = sys.stderr
    bytes_to_sample = 20000
    program_version = {} # keep track of what software we run and the version

    def __init__(self, file_names, platform=None, work_dir=None, interleaved=False):
        """
        create a read_set object
        """
        ReadLibrary.LOG.write("ReadLibrary( %s, platform=%s, interleaved=%s\n"%(file_names, str(platform), str(interleaved)))

        self.num_reads = 0
        self.problem = []
        self.layout = 'na'
        self.format = 'fastq'
        self.platform = 'na'
        self.length_class = 'na'
        self.versions = []
        if platform:
            self.platform = platform
            if platform in ("illumina", "iontorrent"):
                self.length_class = 'short'
            elif platform in ("pacbio", "nanopore"):
                self.length_class = 'long'

        self.transformation ='original'
        self.versions = []
        self.files = []
        input_files = []
        ReadLibrary.LOG.write("read files passed to constructor: {}, type={}\n".format(file_names, type(file_names)))
        if type(file_names) is str:
            input_files.append(file_names)
            ReadLibrary.LOG.write("single-end\n")
        else: # an array of strings
            input_files.extend(file_names)
            ReadLibrary.LOG.write("paired-end\n")
        self.files = [0]*len(input_files)
        self.file_size = [0]*len(input_files)
        for i, read_file in enumerate(input_files):
            if os.path.exists(read_file):
                self.file_size[i] = os.path.getsize(read_file)
                if work_dir: # symlink files to where work will be performed
                    dir_name, file_base = os.path.split(read_file)
                    self.files[i] = file_base
                    ReadLibrary.LOG.write("symlinking %s to %s\n"%(os.path.abspath(read_file), os.path.join(work_dir, file_base)))
                    if os.path.exists(os.path.join(work_dir,file_base)):
                        ReadLibrary.LOG.write("first deleting file {}\n".format(os.path.join(work_dir,file_base)))
                        os.remove(os.path.join(work_dir,file_base))
                    os.symlink(os.path.abspath(read_file), os.path.join(work_dir,file_base))
            else:
                comment = "file does not exist: %s\n"%read_file
                ReadLibrary.LOG.write(comment)
                comment = "self.files = {}\n".format(self.files)
                ReadLibrary.LOG.write(comment)
                self.problem.append(comment)
                raise Exception(comment)

    def store_current_version(self):
        print("store_current_version: self={}, s.v={}".format(self, self.versions))
        current_version = copy.deepcopy(self)
        self.versions.append(current_version)
        print("cv={}, s.v={}".format(current_version, self.versions))

    def bunzip_reads(self):
        ReadLibrary.LOG.write("bunzip_reads()\n")
        self.store_current_version()
        self.transformation="bunzip2"
        startTime = time()
        for i, read_file in enumerate(self.files):
            if read_file.endswith('.bz2'):
                uncompressed_file = read_file[:-4] # trim off '.bz2''
                uncompressed_file = os.path.basename(uncompressed_file) # will write to current working directory
                self.files[i] = uncompressed_file
                with open(os.path.join(uncompressed_file), 'wb') as OUT:
                    with open(os.path.join(read_file), 'rb') as IN:
                        OUT.write(bz2.decompress(IN.read()))
                comment = "decompressing bz2 file %s to %s"%(read_file, uncompressed_file)
                ReadLibrary.LOG.write(comment+"\n")
                self.file_size[i] = os.path.getsize(uncompressed_file)
            else:
                comment = "file {} does not end in '.bz2', not decompressing."
                ReadLibrary.LOG.write(comment)
                self.transformation = comment

        self.command = "python bz2.decompress()"
        self.processing_time = time() - startTime
        ReadLibrary.LOG.write("bunzip_reads duration: {}\n".format(self.processing_time))
        return

    def trim_short_reads(self):
        startTime = time()
        ReadLibrary.LOG.write("trim_short_reads()\n")

        command = ["trim_galore", "--version"]
        proc = subprocess.run(command, shell=False, capture_output=True, text=True)
        version_text = proc.stdout.strip()
        m = re.search(r"(version\s+\S+)", version_text)
        if m:
            trim_galore_version = "trim_galore " + m.group(1)
            ReadLibrary.program_version['trim_galore'] = trim_galore_version

        self.problem = []
        read_file_base = re.sub("(.*?)\..*", "\\1", self.files[0])
        read_file_base = re.sub("(.*)_R[12].*", "\\1", read_file_base)
        trim_directory = read_file_base + "_trim_dir"
        command = ['trim_galore', '-j', str(ReadLibrary.NUM_THREADS), '-o', trim_directory]
        if len(self.files) > 1:
            command.extend(["--paired", self.files[0], self.files[1]])
        else:
            command.append(self.files[0])

        ReadLibrary.LOG.write("command: "+" ".join(command)+"\n")
        self.command = command
        proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE, text=True)
        trimGaloreStderr = proc.stderr.read()
        return_code = proc.wait()
        ReadLibrary.LOG.write("return code = %d\n"%return_code)
        trimReads = glob.glob(trim_directory + "/*val_?.fq")
        if not trimReads:
            trimReads = glob.glob(trim_directory + "/*fq")
        if not trimReads:
            trimReads = glob.glob(trim_directory + "/*fq.gz")
        if trimReads:
            print("trimmed reads = "+str(trimReads))
            self.store_current_version()
            for i, read_file in enumerate(trimReads):
                new_read_file = os.path.basename(read_file)
                shutil.move(read_file, new_read_file)
                self.file_size[i] = os.path.getsize(new_read_file)
                self.files[i] = new_read_file

            self.study_reads()
            self.transformation = "trim-galore"

            report_files = glob.glob(trim_directory +'/*trimming_report.txt')
            if report_files:
                new_report_file = read_file_base + "_trimming_report.txt"
                with open(new_report_file, 'w') as F:
                    for rf in report_files:
                        with open(rf) as RF:
                            F.write(RF.read())
                self.trim_report = new_report_file
        else:
            comment = "trim_galore did not generate trimmed reads"
            ReadLibrary.LOG.write(comment+"\n")
            self.problem.append(comment)


        self.processing_time = time() - startTime
        ReadLibrary.LOG.write("trim_short_reads duration: {}\n".format(self.processing_time))
        return

    def saveTrimReport(self, save_dir):
        pass
        #for name in sorted(self.read_set):
        #    for read_version in self.read_set[name]['versions']:
        #        if 'trim report' in read_version:
        #            shutil.copy(read_version['trim report'], save_dir)

    def study_reads(self):
        """
        Determine read count and avg read length. 
        If paired, read both files.  
        """
        startTime = time()
        ReadLibrary.LOG.write("\nstart study_reads()\n")
        # see if we need to handle bz2 compression - some programs cannot handle it
        if self.files[0].endswith(".bz2"):
            self.bunzip_reads()
        
        self.format = 'fastq'
        self.avg_length = 0
        self.avg_quality = 0
        self.num_reads = 0
        self.num_bases = 0
        self.problem = []
        ReadLibrary.LOG.write("file(s): "+':'.join(self.files)+"\n")

        file1 = self.files[0]
        file2 = None
        if not os.path.exists(file1):
            print("file {} does not exist".format(file1))
            print("cur dir = {}".format(os.getcwd()))
            print("dir listing: {}".format(os.listdir('.')))
            raise Exception("file {} does not exist".format(file1))
        self.file_size[0] = os.path.getsize(file1)
        if len(self.files) > 1:
            self.layout = 'paired-end'
            file2 = self.files[1]
            self.file_size[1] = os.path.getsize(file2)
        else:
            self.layout = 'single-end'
        if file1.endswith("gz"):
            F1 = gzip.open(file1, 'rt')
            if file2:
                F2 = gzip.open(file2, 'rt')
        elif file1.endswith("bz2"):
            F1 = bz2.BZ2File(file1, 'rt')
            if file2:
                F2 = bz2.BZ2File(file2, 'rt')
        else:
            F1 = open(file1, 'rt')
            if file2:
                F2 = open(file2, 'rt')

        line = str(F1.readline().rstrip())
        ReadLibrary.LOG.write("in study_reads, first line of {} is {}\n".format(file1, line))
        sample_read_id = line.split(' ')[0]
        F1.seek(0)

        read_ids_paired = True
        totalReadLength = 0
        maxReadLength = 0
        sumQuality = 0
        numQualityPositionsSampled = 0
        numQualityLinesSampled = 0
        qualityPositionsToSamplePerRead = 50
        numReadsToSampleForQuality = 10000
        readNumber = 0
        seqLen = 0

        if sample_read_id.startswith('>'):
            ReadLibrary.LOG.write("read ID starts with >: "+sample_read_id)
            self.format = 'fasta'

            for line in F1:
                if line.startswith(">"):
                    readNumber += 1
                    if seqLen > maxReadLength:
                        maxReadLength = seqLen
                    totalReadLength += seqLen
                    seqLen = 0
                else:
                    seqLen += len(line)-1

        else: # fastq
            for i, line1 in enumerate(F1):
                if file2:
                    line2 = F2.readline()
                    if not line2:
                        comment = "Number of reads differs between {} and {} after {}".format(file1, file2, readNumber)
                        ReadLibrary.LOG.write(comment+"\n")
                        self.problem.append(comment)
                        break
                    if False and i % 4 == 0 and read_ids_paired:
                        read_id_1 = str(line1).split(' ')[0] # get part up to first space, if any 
                        read_id_2 = str(line2).split(' ')[0] # get part up to first space, if any 
                        if not read_id_1 == read_id_2:
                            diff = findSingleDifference(read_id_1, read_id_2)
                            if diff == None or sorted((read_id_1[diff[0]:diff[1]], read_id_2[diff[0]:diff[1]])) != ('1', '2'):
                                read_ids_paired = False
                                self.problem.append("id_mismatch at read %d: %s vs %s"%(readNumber+1, read_id_1, read_id_2))
                if i % 4 == 1:
                    seqLen = len(line1)-1
                    totalReadLength += seqLen
                    maxReadLength = max(maxReadLength, seqLen) 
                    readNumber += 1
                    if file2:
                        seqLen2 = len(line2)-1
                        totalReadLength += seqLen2
                        maxReadLength = max(maxReadLength, seqLen2) 
                        readNumber += 1
                if i % 4 == 3 and numQualityLinesSampled < numReadsToSampleForQuality:
                    numQualityLinesSampled += 1
                    j = 0
                    for j, qual in enumerate(line1[:qualityPositionsToSamplePerRead]):
                        sumQuality += ord(qual) - 33
                        numQualityPositionsSampled += 1
                if False and readNumber % 100000 == 0:
                    ReadLibrary.LOG.write("number of reads and bases studied so far: \t{}\t{}\n".format(readNumber, totalReadLength))
                    ReadLibrary.LOG.flush()

        F1.close()
        if file2:
            F2.close()

        avgReadLength = 0
        if readNumber:
            avgReadLength = totalReadLength/readNumber
        if file2:
            avgReadLength/=2
        avgReadQuality = 0
        if numQualityPositionsSampled:
            avgReadQuality = sumQuality / float(numQualityPositionsSampled)
        ReadLibrary.LOG.write("avgReadLength={}, avgReadQuality={}, maxLength={}, numReads={}, numBases={}\n".format(avgReadLength, avgReadQuality, maxReadLength, readNumber, totalReadLength))
        self.avg_length = avgReadLength
        self.max_read_len = maxReadLength
        self.num_reads = readNumber
        self.num_bases = totalReadLength
        self.sample_read_id = sample_read_id 
        self.avg_quality = avgReadQuality 

        if not self.platform in ('illumina', 'iontorrent', 'nanopore', 'pacbio'):
            ReadLibrary.LOG.write("platform = {}, need to inferPlatform\n".format(self.platform))
            platform = inferPlatform(sample_read_id, maxReadLength, avgReadQuality)
            self.platform = platform
            ReadLibrary.LOG.write("platform inferred to be {}\n".format(platform))

        self.length_class = ["short", "long"][maxReadLength >= ReadLibrary.MAX_SHORT_READ_LENGTH]
        if file2 and maxReadLength >= ReadLibrary.MAX_SHORT_READ_LENGTH:
            comment = "paired reads appear to be long, expected short: {}".format(",".join(self.files))
            ReadLibrary.LOG.write(comment+"\n")
            self.problem.append(comment)
        ReadLibrary.LOG.write("analysis of {} shows read_number = {} and total_bases = {}\n".format(":".join(self.files), readNumber, totalReadLength))
        ReadLibrary.LOG.write("duration of study_reads was %d seconds\n"%(time() - startTime))

        return

    def normalize_read_depth(self, target_depth=100):
        #use BBNorm
        startTime = time()
        comment = "normalize read depth using BBNorm"
        ReadLibrary.LOG.write(comment+"\n")

        file1 = self.files[0]
        out_file1 = os.path.basename(file1) # will write to current working directory
        file2 = None
        if len(self.files) > 1:
            file2 = self.files[1]
            out_file2 = os.path.basename(file2) # will write to current working directory
        suffix = "_normalized.fastq"

        if file1.endswith("gz"):
            out_file1 = file1[:-3]
            if file2 and file2.endswith("gz"):
                out_file2 = file2[:-3]
        if file1.endswith(".fq"):
            out_file1 = file1[:-3]
            if file2 and file2.endswith(".fq"):
                out_file2 = file2[:-3]
        if file1.endswith(".fastq"):
            out_file1 = file1[:-6]
            if file2 and file2.endswith(".fastq"):
                out_file2 = file2[:-6]
        bbnorm_stdout = out_file1 + "_bbnorm_stats.txt"
        out_file1 = out_file1 + suffix
        if file2:
            out_file2 = out_file2 + suffix

        bbnorm_fh = open(bbnorm_stdout, 'w')
        command = ['bbnorm.sh', 'in='+file1, 'out='+out_file1]
        if file2:
            command.extend(('in2='+file2, 'out2='+out_file2))
        command.extend(['threads={}'.format(ReadLibrary.NUM_THREADS), '-Xmx{:.0f}g'.format(ReadLibrary.MEMORY * 0.85)])

        ReadLibrary.LOG.write("normalize, command line = "+" ".join(command)+"\n")
        self.command = " ".join(command)
        proc = subprocess.Popen(command, shell=False, stderr=bbnorm_fh)
        proc.wait()
        bbnorm_fh.close()

        if os.path.exists(out_file1) and os.path.getsize(out_file1) > 0:
            # if successful, store current version and move on to normalized version
            self.store_current_version()
            self.files[0]=out_file1 
            if file2:
                self.files[1] = out_file2
            self.study_reads()    
            comment = "normalize read depth using BBNorm"
            self.transformation=comment
            self.process_time = time() - startTime

            proc = subprocess.run("bbnorm.sh", capture_output=True, text=True)
            for line in proc.stdout.split("\n"):
                m = re.match("Last modified\s*(\S.*\S)", line)
                if m:
                    bbnorm_version = m.group(1)
                    ReadLibrary.program_version['bbnorm'] = bbnorm_version
                    break

        else:
            # keep current, unnormalized, version
            ReadLibrary.LOG.write("bbnorm failed to normalize\n")
            

        ReadLibrary.LOG.write("normalize process time: {}\n".format(time() - startTime))
        return

    def down_sample_reads(self, max_bases=0):
        """
        read file over size limit, down-sample using seqtk
        """
        startTime = time()
        ReadLibrary.LOG.write("down_sample_reads()\n")
        self.store_current_version()

        if not max_bases:
            max_bases = MAX_BASES
        prop_to_sample = float(max_bases) / self.num_bases
        if prop_to_sample > 1:
            ReadLibrary.LOG.write("down_sample_reads calculated prop_to_sample to be over 1 ({})".format(prop_to_sample))
            return
        comment = "down-sample by {:.3f}X to approximately {} bases".format(prop_to_sample, max_bases)
        self.transformation = comment
        self.command = ''
        ReadLibrary.LOG.write(comment+"\n")
        ReadLibrary.LOG.write("files = "+", ".join(self.files)+"\n")
        
        suffix = "_sampled.fq"
        for i, read_file in enumerate(self.files):
            out_file = os.path.basename(read_file) # will write to current working directory
            if out_file.endswith("gz"):
                out_file = out_file[:-3]
            if out_file.endswith(".fq"):
                out_file = out_file[:-3]
            if out_file.endswith(".fastq"):
                out_file = out_file[:-6]
            out_file += suffix

            self.files[i] = out_file
            
            out_fh = open(out_file, 'w')
            command = ['seqtk', 'sample', read_file, '{:.3f}'.format(prop_to_sample)] 
            ReadLibrary.LOG.write("downsample, command line = "+" ".join(command)+"\n")
            self.command = " ".join(command)+"\n"
            proc = subprocess.Popen(command, shell=False, stdout=out_fh)
            proc.wait()
            out_fh.close()
            self.file_size[i] = os.path.getsize(out_file)

        proc = subprocess.run("seqtk", capture_output=True, text=True)
        for line in proc.stderr.split("\n"):
            m = re.match("Version:\s*(\S.*\S)", line)
            if m:
                seqtk_version = m.group(1)
                ReadLibrary.program_version['seqtk'] = seqtk_version
                break

        ReadLibrary.LOG.write("after: files = "+", ".join(self.files)+"\n")
        self.process_time = time() - startTime
        self.study_reads()    
        ReadLibrary.LOG.write("duration of down_sample_reads and study_reads: %d seconds\n"%(time() - startTime))
        return
    
    def writeHtmlSection(self, HTML):
        versions = list(self.versions)
        versions.append(self)
        print("versions = {}".format(versions))
        print("self.versions = {}".format(self.versions))
        #versions.reverse()
        print("versions = {}".format(versions))
        #HTML.write("<p>"+name+"</p><br>\n")
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th>#</th><th>Process</th><th>File(s)</th><th>Num Reads</th><th>Million Bases</th><th>Avg Length</th></tr></thead>
            <tbody>
            """)
        for i, read_version in enumerate(versions):
            print("html for {}".format(read_version))
            HTML.write("<tr><td>{}</td>".format(i))
            HTML.write("<td>{}</td>".format(read_version.transformation))
            HTML.write("<td>{}</td>".format(" ".join(read_version.files)))
            HTML.write("<td>{}</td>".format(read_version.num_reads))
            HTML.write("<td>{:.2f}</td>".format(read_version.num_bases / 1e6))
            HTML.write("<td>{:.1f}</td>".format(read_version.avg_length))
            HTML.write("</tr>\n")
        HTML.write("</tbody></table>\n")
