#!/usr/bin/env python
import sys
import subprocess
import argparse
import gzip
import bz2
import os
import os.path
import re
import shutil
from time import time, localtime, strftime, sleep

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

def inferPlatform(read_id, maxReadLength, avgReadQuality):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    if read_id.startswith(">"):
        return "fasta"
    if maxReadLength < FastqPreprocessor.MAX_SHORT_READ_LENGTH:
        # example illumina read id
        #@D00553R:173:HG53VBCXY:2:1101:1235:2074 1:N:0:ACAGTGAT
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

class FastqPreprocessor:
    MAX_BASES=5e9
    bytes_to_sample = 20000
    MAX_SHORT_READ_LENGTH = 600
    DEFAULT_THREADS = 4
    DEFAULT_MEMORY = 64 # in gigabytes

    def __init__(self, working_directory=".", log=None, max_bases = 0, threads = None, memory = None):
        self.LOG = sys.stderr # create a log file at start of main()
        self.WORK_DIR = '.'
        if working_directory:
            self.WORK_DIR = working_directory
        if log:
            self.LOG = log
        else:
            self.LOG = open(os.path.join(self.WORK_DIR, "fastq_processor_log.txt"), 'w')
        self.LOG.write("starting FastqPreprocessor\n")
        self.LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
        self.LOG.write("Work directory is "+self.WORK_DIR+"\n\n")
        self.read_set = {}
        self.preprocess_steps = []
        self.program_version = {}
        self.problem = []
        self.max_bases = FastqPreprocessor.MAX_BASES
        if max_bases:
            self.max_bases=max_bases
        self.threads = FastqPreprocessor.DEFAULT_THREADS
        if threads:
            self.threads = threads
        self.memory = FastqPreprocessor.DEFAULT_MEMORY
        if memory:
            self.memory = memory

    def setMaxBases(self, max_bases):
        self.max_bases = max_bases

    def setThreads(self, threads):
        self.threads = threads

    def addProcessStep(self, step):
        self.preprocess_steps.append(step)

    def getPlatformReads(self, platform):
        retval = []
        for name in self.read_set:
            read_set = self.read_set[name]
            if read_set['platform'] == platform:
                retval.append(read_set['versions'][-1]['files'])
        return retval

    def getPairedReads(self):
        retval = []
        for name in self.read_set:
            if self.read_set[name]['length_class'] == 'short':
                files = self.read_set[name]['versions'][-1]['files']
                if len(files) == 2:
                    retval.append(files)
        return retval

    def getUnpairedShortReads(self):
        retval = []
        for name in self.read_set:
            if self.read_set[name]['length_class'] == 'short':
                files = self.read_set[name]['versions'][-1]['files']
                if len(files) == 1:
                    retval.append(files[0])
        return retval

    def getShortReads(self):
        retval = []
        self.LOG.write("getShortReads\n")
        for name in self.read_set:
            self.LOG.write("read_set {}, length_class: {}\n".format(name, self.read_set[name]['length_class']))
            self.LOG.write("read_set explicit: {}\n\n".format(self.read_set[name]))
            if self.read_set[name]['length_class'] == 'short':
                retval.append(self.read_set[name]['versions'][-1]['files'])
                self.LOG.write("     adding files "+", ".join(self.read_set[name]['versions'][-1]['files'])+"\n")
        self.LOG.flush()
        return retval

    def getLongReads(self):
        retval = []
        for name in self.read_set:
            if self.read_set[name]['length_class'] == 'long':
                retval.append(self.read_set[name]['versions'][-1]['files'][0])
        return retval

    def getPacbioReads(self):
        retval = []
        for name in self.read_set:
            if self.read_set[name]['platform'] == 'pacbio':
                retval.append(self.read_set[name]['versions'][-1]['files'][0])
        return retval

    def getNanoporeReads(self):
        retval = []
        for name in self.read_set:
            if self.read_set[name]['platform'] == 'nanopore':
                retval.append(self.read_set[name]['versions'][-1]['files'][0])
        return retval

    def getReadPlatform(self, reads):
        if reads in self.read_set:
            return self.read_set[reads]['platform']
        else:
            return None

    def sampleReads(self, filename):
        self.LOG.write("sampleReads()\n")
        # figures out Read_file_type
        #return read_format and sample of read ids
        read_format = 'na'
        read_id_sample = []

        if filename.endswith("gz"):
            F = gzip.open(filename)
        elif filename.endswith("bz2"):
            F = bz2.BZ2File(filename)
        else:
            F = open(filename)
        text = F.read(FastqPreprocessor.bytes_to_sample) #read X number of bytes for text sample
        F.close()

        self.LOG.write("  file text sample %s:\n%s\n\n"%(filename, text[0:50]))
        lines = text.split("\n")
        readLengths = []
        if len(lines) < 2:
            comment = "in sampleReads for %s: text sample (length %d) lacks at least 2 lines"%(filename, len(text))
            self.LOG.write(comment+"\n")
            self.problem.append(comment)
        if lines[0].startswith("@"):
            read_format = 'fastq'
            for i, line in enumerate(lines):
                if i % 4 == 0:
                    read_id_sample.append(line.split(' ')[0]) # get part up to first space, if any 
                elif i % 4 == 1:
                    readLengths.append(len(line)-1)
        elif lines[0].startswith(">"):
            read_format = 'fasta'
            read_id_sample.append(lines[0].split()[0])
            seq = ""
            for line in lines:
                if line.startswith(">"):
                    readLengths.append(len(seq))
                    seq = ""
                else:
                    seq += line.rstrip()
        max_read_length = 0
        if len(readLengths) > 0:
            max_read_length = max(readLengths)
        else:
            comment = "in sampleReads for %s: text sample (length %d) did not contain any sequences to determine length from."%(filename, len(text))
            self.LOG.write(comment+"\n")
            self.problem.append(comment)
        if len(read_id_sample) > 1:
            read_id_sample = read_id_sample[:-1] # last entry might be truncated, avoid it
        self.LOG.write("read type %s, maximum read length %.1f\n"%(read_format, max_read_length))
        return read_id_sample, max_read_length

    def registerReads(self, reads, platform=None, interleaved=False):
        """
        create a read_set entry in details for these reads
        move read files to working directory to allow relative paths
        """
        self.LOG.write("registerReads( %s, platform=%s, interleaved=%s\n"%(reads, str(platform), str(interleaved)))

        name = reads.split(":")[0]
        name = name.split("/")[-1]
        name = name.split(".")
        while name[-1] in ("fastq", "fq", "gz", "bz2"):
            name = name[:-1]
        name = ".".join(name)


        if name in self.read_set:
            comment = "duplicate registration of reads %s"%name
            self.LOG.write(comment+"\n")
            reads += "_V"
            v_num = 2
            while name+str(v_num) in self.read_set:
                v_num += 1
            name += str(v_num)

        read_set = {} # high-level representation of read set, with specific versions included within (e.g., trimmed version)
        read_set['name'] = name
        self.read_set[name] = read_set

        read_set["problem"] = []
        read_set["layout"] = 'na'
        read_set["platform"] = None
        read_set['length_class'] = 'na'
        read_set['delim'] = ''
        read_set['original_name'] = reads
        read_set['versions'] = []
        read_set['failed'] = []
        if platform:
            read_set['platform'] = platform
            if platform in ("illumina", "iontorrent"):
                read_set['length_class'] = 'short'
            elif platform in ("pacbio", "nanopore"):
                read_set['length_class'] = 'long'

        read_version = {}
        read_version['read_set'] = read_set
        read_version['transformation'] ='original'
        read_version['files'] = []
        read_version['file_size'] = 0
        read_set['versions'].append(read_version)
        delim = ["%", ":"][":" in reads] # ":" if it is, else "%"
        read_list = reads.split(delim)
        if len(read_list) > 1:
            read_set['delim'] = delim
        for read_file in read_list:
            if not os.path.exists(read_file):
                comment = "file does not exist: %s"%read_file
                self.LOG.write(comment+"\n")
                read_set["problem"].append(comment)
                return None
            read_version['file_size'] += os.path.getsize(read_file)
            dir_name, file_base = os.path.split(read_file)
            if os.path.abspath(dir_name) != os.path.abspath(self.WORK_DIR):
                self.LOG.write("symlinking %s to %s\n"%(os.path.abspath(read_file), os.path.join(self.WORK_DIR,file_base)))
                if os.path.exists(os.path.join(self.WORK_DIR,file_base)):
                    self.LOG.write("first deleting file {}\n".format(os.path.join(self.WORK_DIR,file_base)))
                    os.remove(os.path.join(self.WORK_DIR,file_base))
                os.symlink(os.path.abspath(read_file), os.path.join(self.WORK_DIR,file_base))
            read_version['files'].append(file_base)
        return read_version

    def preprocess_reads(self):
        self.LOG.write("preprocess_reads()\n")
        startTime = time()
        prev_dir = os.getcwd()
        os.chdir(self.WORK_DIR)
        for name in sorted(self.read_set):
            self.LOG.write("process read_set {}\n".format(name))
            read_set = self.read_set[name]
            read_version = read_set['versions'][-1]

            # first see if we need to handle bz2 compression - most program can handle it
            if read_version['files'][0].endswith(".bz2"):
                new_read_version = self.bunzip_reads(read_version)
                if new_read_version['file_size'] > 0:
                    read_set['versions'].append(new_read_version)
                    read_version = new_read_version
                else:
                    self.LOG.write("bunzip_reads yielded no output")
                    read_set['failed'].append(new_read_version)

            self.study_reads(read_version)

            for step in self.preprocess_steps:
                self.LOG.write("\nstep {}\n".format(step))
                if step == 'trim_short_reads' and read_set['length_class'] == 'short':
                    new_read_version = self.trim_short_reads(read_version)
                    if new_read_version['file_size'] > 0:
                        read_set['versions'].append(new_read_version)
                        read_version = new_read_version
                        self.study_reads(read_version)
                    else:
                        self.LOG.write("trim_short_reads yielded no output reads file")
                        read_set['failed'].append(new_read_version)

                elif step == 'normalize':
                    new_read_version = self.normalize_read_depth(read_version)
                    if new_read_version['file_size'] > 0:
                        read_set['versions'].append(new_read_version)
                        read_version = new_read_version
                    else:
                        self.LOG.write("trim_short_reads yielded no output reads file")
                        read_set['failed'].append(new_read_version)

            if self.max_bases and read_version["num_bases"] > self.max_bases:
                comment = "need to downsample reads {}".format(read_set['delim'].join(read_version['files']))
                self.LOG.write(comment+"\n")
                read_version = self.down_sample_reads(read_version, max_bases=self.max_bases)
                read_set['versions'].append(read_version)

        self.LOG.write("preprocessing time in seconds: {}\n".format(time() - startTime))
        os.chdir(prev_dir)
        return 

    def summarize_preprocess(self):
        for name in sorted(self.read_set):
            print("\nread set name: {}\n".format(name))
            read_set = self.read_set[name]
            print("\nread set: {}\n".format(read_set))
            print("number of versions: {}\n".format(len(read_set['versions'])))
            for i, version in enumerate(read_set['versions']):
                print("version {}".format(i))
                for key in version:
                    if key == 'read_set':
                        continue
                    print("  {}: {}".format(key, version[key]))
        print("done showing read sets\n")


    def bunzip_reads(self, read_version):
        self.LOG.write("bunzip_reads()\n")
        startTime = time()
        new_read_version = {}
        new_read_version['read_set'] =  read_version['read_set']
        new_read_version['files'] = []
        new_read_version['file_size'] = 0
        for read_file in read_version['files']:
            if read_file.endswith('.bz2'):
                uncompressed_file = read_file[:-4] # trim off '.bz2''
                new_read_version['files'].append(uncompressed_file)
                with open(os.path.join(WORK_DIR, uncompressed_file), 'w') as OUT:
                    with open(os.path.join(WORK_DIR, read_file)) as IN:
                        OUT.write(bz2.decompress(IN.read()))
                        comment = "decompressing bz2 file %s to %s"%(read_file, uncompressed_file)
                        self.LOG.write(comment+"\n")
                        new_read_version['file_size'] += os.path.getsize(uncompressed_file)
                        #details["pre-assembly transformation"].append(comment)
            else:
                comment = "file {} does not end in '.bz2', not decompressing."
                self.LOG.write(comment)
                new_read_version['transformation'] = comment

        new_read_version['transformation'] = "remove bz2 compression"
        new_read_version['command'] = "python bz2.decompress()"
        new_read_version['processing_time'] = time() - startTime()
        self.LOG.write("bunzip_reads duration: {}\n".format(new_read_version['processing_time']))
        return new_read_version

    def trim_short_reads(self, read_version):
        startTime = time()
        self.LOG.write("trim_short_reads()\n")

        if "trim_galore" not in self.program_version:
            command = ["trim_galore", "--version"]
            proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE)
            version_text = proc.stdout.read().decode().strip()
            proc.wait()
            m = re.search(r"(version\s+\S+)", version_text)
            if m:
                trim_galore_version = "trim_galore " + m.group(1)
                self.program_version['trim_galore'] = trim_galore_version

        new_read_version = {}
        new_read_version['read_set'] =  read_version['read_set']
        new_read_version['file_size']  = 0
        num_threads = self.threads
        command = ['trim_galore', '-j', str(num_threads), '-o', '.']
        if len(read_version['files']) > 1:
            command.extend(["--paired", read_version['files'][0], read_version['files'][1]])
        else:
            command.append(read_version['files'][0])

        self.LOG.write("command: "+" ".join(command)+"\n")
        new_read_version['command'] = command
        proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE, text=True)
        trimGaloreStderr = proc.stderr.read()
        return_code = proc.wait()
        self.LOG.write("return code = %d\n"%return_code)
        trimReads = re.findall(r"Writing validated paired-end read \d reads to (\S+)", trimGaloreStderr)
        if not trimReads:
            trimReads = re.findall(r"Writing final adapter and quality trimmed output to (\S+)", trimGaloreStderr)
        self.LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
        if not trimReads:
            comment = "trim_galore did not name trimmed reads output files in stderr"
            self.LOG.write(comment+"\n")
            new_read_version['problem'].append(comment)
            return new_read_version
        new_read_version['files']  = []
        for trimmed_read_file in trimReads:
            if re.search("val_[12].fq", trimmed_read_file):
                val_file = trimmed_read_file
                trimmed_read_file = re.sub("val_[12].fq", "trimmed.fastq", val_file)
                shutil.move(val_file, trimmed_read_file)
            if re.search("trimmed.fq", trimmed_read_file):
                temp_file = trimmed_read_file
                trimmed_read_file = trimmed_read_file.replace(".fq", ".fastq")
                shutil.move(temp_file, trimmed_read_file)
            new_read_version['files'].append(trimmed_read_file)
            new_read_version['file_size'] += os.path.getsize(trimmed_read_file)
        comment = "trim_galore, input %s, output %s"%(":".join(read_version['files']), ":".join(trimReads))
        self.LOG.write(comment+"\n")
        new_read_version['transformation'] = comment
        #new_read_version['report'] = 
        m = re.search(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
        if m:
            report_file = m.group(1)
            self.LOG.write("re.search for trim reports returned %s\n"%str(report_file))
            new_read_version["trim report"] = report_file
        new_read_version['processing_time'] = time() - startTime
        self.LOG.write("trim_short_reads duration: {}\n".format(new_read_version['processing_time']))

        return new_read_version

    def study_reads(self, read_version):
        """
        Determine read count and avg read length. Update read_version structure.
        If paired, read both files.  
        """
        startTime = time()
        self.LOG.write("\nstart study_reads()\n")
        read_version['avg_length'] = 0
        read_version['avg_quality'] = 0
        read_version['num_reads'] = 0
        read_version['num_bases'] = 0
        read_version['file_size'] = 0
        self.LOG.write("file(s): "+':'.join(read_version['files'])+"\n")

        file1 = read_version['files'][0]
        file2 = None
        if not os.path.exists(file1):
            print("file {} does not exist".format(file1))
            print("cur dir = {}".format(os.getcwd()))
            print("dir listing: {}".format(os.listdir('.')))
            raise Exception("file {} does not exist".format(file1))
        read_version['file_size'] = os.path.getsize(file1)
        if len(read_version['files']) > 1:
            read_version['layout'] = 'paired-end'
            file2 = read_version['files'][1]
            read_version['file_size'] += os.path.getsize(file2)
        else:
            read_version['layout'] = 'single-end'
        if file1.endswith("gz"):
            F1 = gzip.open(file1)
            if file2:
                F2 = gzip.open(file2)
        elif file1.endswith("bz2"):
            F1 = bz2.BZ2File(file1)
            if file2:
                F2 = bz2.BZ2File(file2)
        else:
            F1 = open(file1)
            if file2:
                F2 = open(file2)

        line = str(F1.readline().rstrip())
        self.LOG.write("in study_reads, first line of {} is {}\n".format(file1, line))
        sample_read_id = line.split(' ')[0]
        F1.seek(0)
        if sample_read_id.startswith('>'):
            raise Exception("read ID starts with >: "+sample_read_id)

        read_ids_paired = True
        totalReadLength = 0
        maxReadLength = 0
        sumQuality = 0
        numQualityPositionsSampled = 0
        qualityPositionsToSamplePerRead = 150
        numReadsToSampleForQuality = 10000

        readNumber = 0
        for i, line1 in enumerate(F1):
            if file2:
                line2 = F2.readline()
                if not line2:
                    comment = "Number of reads differs between {} and {} after {}".format(file1, file2, readNumber)
                    self.LOG.write(comment+"\n")
                    read_version["problem"].append(comment)
                    break
                if False and i % 4 == 0 and read_ids_paired:
                    read_id_1 = str(line1).split(' ')[0] # get part up to first space, if any 
                    read_id_2 = str(line2).split(' ')[0] # get part up to first space, if any 
                    if not read_id_1 == read_id_2:
                        diff = findSingleDifference(read_id_1, read_id_2)
                        if diff == None or sorted((read_id_1[diff[0]:diff[1]], read_id_2[diff[0]:diff[1]])) != ('1', '2'):
                            read_ids_paired = False
                            read_version["problem"].append("id_mismatch at read %d: %s vs %s"%(readNumber+1, read_id_1, read_id_2))
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
            if i % 4 == 3 and readNumber < numReadsToSampleForQuality:
                for j, qual in enumerate(line1[:qualityPositionsToSamplePerRead]):
                    sumQuality += ord(qual) - 33
                    numQualityPositionsSampled += 1
            if False and readNumber % 100000 == 0:
                self.LOG.write("number of reads and bases studied so far: \t{}\t{}\n".format(readNumber, totalReadLength))
                self.LOG.flush()

        F1.close()
        if file2:
            F2.close()

        avgReadLength = 0
        if readNumber:
            avgReadLength = totalReadLength/readNumber
        if file2:
            avgReadLength/=2
        avgReadQuality = sumQuality / float(numQualityPositionsSampled)
        self.LOG.write("avgReadLength={}, avgReadQuality={}, maxLength={}, numReads={}, numBases={}\n".format(avgReadLength, avgReadQuality, maxReadLength, readNumber, totalReadLength))
        read_version['avg_length'] = avgReadLength
        read_version['max_read_len'] = maxReadLength
        read_version['num_reads'] = readNumber
        read_version['num_bases'] = totalReadLength
        read_version['sample_read_id'] = sample_read_id 
        read_version['avg_quality'] = avgReadQuality 

        if not read_version['read_set']['platform'] in ('illumina', 'iontorrent', 'nanopore', 'pacbio'):
            self.LOG.write("platform = {}, need to inferPlatform\n".format(read_version['read_set']['platform']))
            platform = inferPlatform(sample_read_id, maxReadLength, avtReadQuality)
            read_version['read_set']['platform'] = platform
            self.LOG.write("platform inferred to be {}\n".format(platform))

        read_version['read_set']['length_class'] = ["short", "long"][maxReadLength >= FastqPreprocessor.MAX_SHORT_READ_LENGTH]
        if file2 and maxReadLength >= FastqPreprocessor.MAX_SHORT_READ_LENGTH:
            comment = "paired reads appear to be long, expected short: %s"%read_version['delim'].join(read_version['files'])
            self.LOG.write(comment+"\n")
            read_version['problem'].append(comment)
        self.LOG.write("analysis of {} shows read_number = {} and total_bases = {}\n".format(":".join(read_version['files']), readNumber, totalReadLength))
        self.LOG.write("duration of study_reads was %d seconds\n"%(time() - startTime))

        read_version['study_reads_time'] = time() - startTime
        return read_version

    def normalize_read_depth(self, read_version, target_depth=100):
        #use BBNorm
        startTime = time()
        self.LOG.write("normalize_read_depth()\n")
        comment = "normalize read depth using BBNorm"
        self.LOG.write(comment+"\n")
        new_read_version = {}
        new_read_version['transformation'] = comment

        new_read_version['read_set'] = read_version['read_set']
        file1 = read_version['files'][0]
        file2 = None
        if len(read_version['files']) > 1:
            file2 = read_version['files'][1]
        suffix = "_normalized.fastq"

        out_file1, out_file2 = file1, file2
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

        new_read_version['files']=[out_file1]
        if file2:
            new_read_version['files'].append(out_file2)
        
        bbnorm_fh = open(bbnorm_stdout, 'w')
        command = ['bbnorm', '-Xmx{}g'.format(self.memory), 'threads={}'.format(self.threads), 'in='+file1, 'out='+out_file1]
        if file2:
            command.extend(('in2='+file2, 'out2='+out_file2))
        command.append("threads={}".format(self.threads))

        self.LOG.write("normalize, command line = "+" ".join(command)+"\n")
        new_read_version['command'] = " ".join(command)
        proc = subprocess.Popen(command, shell=False, stderr=bbnorm_fh)
        proc.wait()
        bbnorm_fh.close()
        new_read_version['process_time'] = time() - startTime
        self.LOG.write("normalize process time: {}\n".format(new_read_version['process_time']))

        self.study_reads(new_read_version)    

        return new_read_version

    def down_sample_reads(self, read_version, max_bases=0):
        """
        read file over size limit, down-sample using seqtk
        """
        startTime = time()
        self.LOG.write("down_sample_reads()\n")

        if not max_bases:
            max_bases = FastqPreprocessor.MAX_BASES
        prop_to_sample = float(max_bases) / read_version['num_bases']
        self.LOG.write("num_bases = {}, max_bases= {}, prop to sample = {:.3f}\n".format(read_version['num_bases'], max_bases, prop_to_sample))
        comment = "down-sample by {:.3f} to approximately {} bases".format(prop_to_sample, self.max_bases)
        new_read_version = {}
        new_read_version['transformation'] = comment
        new_read_version['file_size'] = 0
        new_read_version['files'] = []
        new_read_version['command'] = ''
        self.LOG.write(comment+"\n")
        
        if prop_to_sample > 1:
            self.LOG.write("down_sample_reads calculated prop_to_sample to be over 1 ({})".format(prop_to_sample))
            #raise Exception("down_sample_reads calculated prop_to_sample to be over 1 ({})".format(prop_to_sample))
            return new_read_version
        new_read_version['read_set'] = read_version['read_set']
        suffix = "_downsampled.fastq"
        for read_file in read_version['files']:
            out_file = read_file
            if read_file.endswith("gz"):
                out_file = read_file[:-3]
            if read_file.endswith(".fq"):
                out_file = read_file[:-3]
            if read_file.endswith(".fastq"):
                out_file = read_file[:-6]
            out_file = out_file + suffix

            new_read_version['files'].append(out_file)
            
            out_fh = open(out_file, 'w')
            command = ['seqtk', 'sample', read_file, '{:.3f}'.format(prop_to_sample)] 
            self.LOG.write("downsample, command line = "+" ".join(command)+"\n")
            new_read_version['command'] += " ".join(command)+"\n"
            proc = subprocess.Popen(command, shell=False, stdout=out_fh)
            proc.wait()
            out_fh.close()
            new_read_version['file_size'] += os.path.getsize(out_file)

        new_read_version['process_time'] = time() - startTime
        self.study_reads(new_read_version)    
        self.LOG.write("duration of down_sample_reads and study_reads: %d seconds\n"%(time() - startTime))
        return new_read_version
         
    def writeHtmlSection(self, HTML):
        HTML.write('<div>\n<h2>Fastq Preprocessing</h2>\n')
        for name in sorted(self.read_set):
            HTML.write("<p>"+name+"</p><br>\n")
            HTML.write("""
            <table class="med-table kv-table">
                <thead class="table-header">
                <tr> <th>#</th><th>Process</th><th>File(s)</th><th>Size</th><th>Reads</th><th>Bases</th><th>Time</th></tr></thead>
                <tbody>
                """)
            for i, read_version in enumerate(self.read_set[name]['versions']):
                HTML.write("<tr><td>{}</td>".format(i))
                HTML.write("<td>{}</td>".format(read_version['transformation']))
                HTML.write("<td>{}</td>".format(",".join(read_version['files'])))
                HTML.write("<td>{}</td>".format(read_version['file_size']))
                HTML.write("<td>{}</td>".format(read_version['num_reads']))
                HTML.write("<td>{:.1f}</td>".format(read_version['avg_length']))
                if 'process_time' in read_version:
                    HTML.write("<td>{}</td>".format(read_version['process_time']))
                if 'study_reads_time' in read_version:
                    HTML.write("<td>{}</td>".format(read_version['study_reads_time']))
                HTML.write("</tr>\n")
        HTML.write("</tbody></table>\n")
        HTML.write("</div\n")

    def writeHtmlReport(self, file_name="fastq_preprocess.html"):
        HTML = open(file_name, 'w')
        HTML.write("<html>\n<H1>Fastq Preprocessing Report</H1>\n")
        self.writeHtmlSection(HTML)
        HTML.write("</html>\n")

if __name__ == "__main__":
    main_return_code = 1 # set to zero when we have an assembly
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='.')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    #parser.add_argument('--sra', metavar='run accessions', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--max_bases', type=int, default=FastqPreprocessor.MAX_BASES, help='downsample reads if more than this total bases.')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=150, help='RAM limit in Gb')
    parser.add_argument('--trim_short_reads', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--normalize', action='store_true', help='normalize read depth with BBNorm at default settings')
    parser.add_argument('--path-prefix', help="Add the given directories to the PATH", nargs='*', required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    WORK_DIR = args.outputDirectory #"p3_assembly" 
    if not os.path.exists(WORK_DIR):
        #shutil.rmtree(WORK_DIR)
        os.mkdir(WORK_DIR)
    WORK_DIR = os.path.abspath(WORK_DIR)
 
    if args.path_prefix:
        os.environ["PATH"] = ":".join(args.path_prefix) + ":" + os.environ["PATH"]

    LOG = open("fastq_processing.log", "w")
    fqProc = FastqPreprocessor(working_directory = args.outputDirectory, log=LOG, memory=args.memory)
    if args.trim_short_reads:
        fqProc.addProcessStep("trim_short_reads")
    if args.normalize:
        fqProc.addProcessStep("normalize")

    if args.anonymous_reads:
        for item in args.anonymous_reads:
            fqProc.registerReads(item, platform = 'anonymous')

    if args.sra:
        for item in args.sra:
            fqProc.registerReads(item, platform = 'sra')

    if args.illumina:
        platform = 'illumina'
        for item in args.illumina:
            interleaved = args.interleaved and item in args.interleaved
            fqProc.registerReads(item, platform=platform, interleaved=interleaved)

    if args.iontorrent:
        platform = 'iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            fqProc.registerReads(item, platform=platform, interleaved=interleaved)

    if args.pacbio:
        for item in args.pacbio:
            fqProc.registerReads(item, platform = 'pacbio')

    if args.nanopore:
        for item in args.nanopore:
            fqProc.registerReads(item, platform = 'nanopore')

    fqProc.preprocess_reads()

    fqProc.summarize_preprocess()

    fqProc.writeHtmlReport("fastq_preprocess_report.html")
