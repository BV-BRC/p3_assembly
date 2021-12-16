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
try:
    import urllib.request as urllib2 # python3
except ImportError:
    import urllib2 #python2
from time import time, localtime, strftime
import json
import glob
#import sra_tools
#import pdb

"""
This script organizes a command line for either 
Unicycler, or canu as appropriate: 
    canu if only long reads (pacbio or nanopore), 
    Unicycler if illumina or iontorrent 
    and Unicycler for hybrid assemblies, eg Illumina plus PacBio
    or Spades if requested
It can auto-detect read types (illumina, iontorrent, pacbio, nanopore)
It can run trim_galore prior to assembling.
It can run Quast to generate assembly quality statistics.
TODO: properly handle different kinds of pacbio reads
TODO: verify that read type identification works in general
"""

DEFAULT_GENOME_SIZE = "5m"
Default_bytes_to_sample = 20000
MAX_SHORT_READ_LENGTH = 600
Read_id_sample = {}
Read_file_type = {}
Avg_read_length = {}
LOG = None # create a log file at start of main()
START_TIME = None
WORK_DIR = None
SAVE_DIR = None
DETAILS_DIR = None

def registerReads(reads, details, platform=None, interleaved=False, supercedes=None, max_bases=None):
    """
    create an entry in details for these reads
    move read files to working directory to allow relative paths
    """
    LOG.write("registerReads( %s, platform=%s, interleaved=%s, supercedes=%s\n"%(reads, str(platform), str(interleaved), str(supercedes)))
    if reads in details["original_items"]:
        comment = "duplicate registration of reads %s"%reads
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    details['original_items'].append(reads)
    
    read_struct = {}
    read_struct['files'] = []
    #read_struct["path"] = []
    read_struct["problem"] = []
    read_struct["layout"] = 'na'
    read_struct["platform"] = 'na'
    read_struct['length_class'] = 'na'
    read_struct['delim'] = ''
    file1, file2 = None, None
    if ":" in reads or "%" in reads:
        if ":" in reads:
            delim = ["%", ":"][":" in reads] # ":" if it is, else "%"
        read1, read2 = reads.split(delim)
        read_struct["delim"] = delim
        if not os.path.exists(read1):
            comment = "file does not exist: %s"%read1
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        if not os.path.exists(read2):
            comment = "file does not exist: %s"%read2
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        dir1, file1 = os.path.split(read1)
        dir2, file2 = os.path.split(read2)
        if os.path.abspath(dir1) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(read1, os.path.join(WORK_DIR,file1)))
            if os.path.exists(os.path.join(WORK_DIR,file1)):
                LOG.write("first deleting file {}\n".format(os.path.join(WORK_DIR,file1)))
                os.remove(os.path.join(WORK_DIR,file1))
            os.symlink(os.path.abspath(read1), os.path.join(WORK_DIR,file1))
        if os.path.abspath(dir2) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(read2, os.path.join(WORK_DIR,file2)))
            if os.path.exists(os.path.join(WORK_DIR,file2)):
                LOG.write("first deleting file {}\n".format(os.path.join(WORK_DIR,file2)))
                os.remove(os.path.join(WORK_DIR,file2))
            os.symlink(os.path.abspath(read2), os.path.join(WORK_DIR,file2))
        read_struct['files'].append(file1)
        read_struct['files'].append(file2)
    else:
        # no ':' or '%' delimiter, so a single file
        if not os.path.exists(reads):
            comment = "file does not exist: %s"%reads
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        if interleaved:
            read_struct["interleaved"] = True
        dir1, file1 = os.path.split(reads)
        if os.path.abspath(dir1) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(reads, os.path.join(WORK_DIR,file1)))
            os.symlink(os.path.abspath(reads), os.path.join(WORK_DIR,file1))
        read_struct['files'].append(file1)
        #read_struct["path"].append(reads)
    read_struct = study_reads(read_struct)
    if False and read_struct["platform"] == 'fasta':
        comment = "sequences are FASTA, cannot process"
        LOG.write(comment)
        sys.stdout.write(comment)
        sys.exit(1)
    too_long = max_bases and read_struct["num_bases"] > max_bases 
    if too_long or ("seq_qual_lengths_differ" in read_struct and read_struct["seq_qual_lengths_differ"]):
        comment = "need to re-write files {}".format(read_struct['delim'].join(read_struct['files']))
        LOG.write(comment+"\n")
        sys.stderr.write(comment+"\n")
        read_struct = rewrite_reads(read_struct, max_bases=max_bases)
        comment = "max_bases limit exceeded, trucating to approximately {} bases".format(max_bases)
        if 'problem' not in read_struct:
            read_struct['problem'] = []
        read_struct['problem'].append(comment)
        LOG.write(comment+"\n")
        #details["pre-assembly transformation"].append(comment)
        file1 = read_struct['files'][0]
        if file2:
            file2 = read_struct['files'][1]
            
    if file1.endswith(".bz2"):
        uncompressed_file1 = file1[:-4]
        with open(os.path.join(WORK_DIR, uncompressed_file1), 'w') as OUT:
            with open(os.path.join(WORK_DIR, file1)) as IN:
                OUT.write(bz2.decompress(IN.read()))
                comment = "decompressing bz2 file %s to %s"%(file1, uncompressed_file1)
                LOG.write(comment+"\n")
                #details["pre-assembly transformation"].append(comment)
                read_struct['files'][0] = uncompressed_file1
    if file2 and file2.endswith(".bz2"):
        uncompressed_file2 = file2[:-4]
        with open(os.path.join(WORK_DIR, uncompressed_file2), 'w') as OUT:
            with open(os.path.join(WORK_DIR, file2)) as IN:
                OUT.write(bz2.decompress(IN.read()))
                comment = "decompressing bz2 file %s to %s"%(file2, uncompressed_file2)
                LOG.write(comment+"\n")
                #details["pre-assembly transformation"].append(comment)
                read_struct['files'][1] = uncompressed_file2

    registeredName = read_struct['delim'].join(read_struct['files'])
    if registeredName in details['reads']:
        comment = "registered name %s already in details[reads]"%registeredName
        LOG.write(comment+"\n")
        #details["problem"].append(comment)

    if supercedes:
        details['derived_reads'].append(registeredName)
        read_struct['supercedes'] = supercedes
        if supercedes in details['reads']:
            details['reads'][supercedes]['superceded_by'] = registeredName
            platform = details['reads'][supercedes]['platform']
            read_struct['platform'] = platform
            try: # swap item in platform[] with this new version
                index = details['platform'][platform].index(supercedes)
                details['platform'][platform][index] = registeredName
            except ValueError:
                comment = "Problem: superceded name %s not found in details_%s"%(registeredName, platform)
    else:
        if platform:
            read_struct['platform'] = platform
            if platform not in details['platform']:
                details["platform"][platform] = []
            details["platform"][platform].append(registeredName)

    details['reads'][registeredName]=read_struct
    return registeredName

def parseJsonParameters(args):
    """ Not fully implemented: should read assembly2 service json file """
    if not os.path.exists(args.params_json):
        raise Exception("cannot find json parameters file %s\n"%args.params_json)
    LOG.write("parseJsonParameters() time=%d\n"%time())
    with open(args.params_json) as json_file:
        data = json.load(json_file)
        #args.output_dir = data['output_file']
        args.min_contig_length = data['min_contig_length']
        if "paired_end_libs" in data:
            for pe in data["paired_end_libs"]:
                pass
        if "single_end_libs" in data:
            for pe in data["single_end_libs"]:
                if pe["platform"] == 'illumina':
                    pass
    return

def study_reads(read_data):
    """
    Determine avg read length. Update details['reads'].
    If paired, read both files. Verify read ID are paired. 
    """
    func_start = time()
    LOG.write("\nstart study_reads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start - START_TIME))
    read_data['avg_len'] = 0
    read_data['length_class'] = 'na'
    read_data['num_reads'] = 0
    read_data['num_bases'] = 0
    LOG.write("file(s): "+read_data['delim'].join(read_data['files'])+"\n")

    file1 = read_data['files'][0]
    file2 = None
    if len(read_data['files']) > 1:
        read_data['layout'] = 'paired-end'
        file2 = read_data['files'][1]
    else:
        read_data['layout'] = 'single-end'
    if file1.endswith("gz"):
        F1 = gzip.open(os.path.join(WORK_DIR, file1))
        if file2:
            F2 = gzip.open(os.path.join(WORK_DIR, file2))
    elif file1.endswith("bz2"):
        F1 = bz2.BZ2File(os.path.join(WORK_DIR, file1))
        if file2:
            F2 = bz2.BZ2File(os.path.join(WORK_DIR, file2))
    else:
        F1 = open(os.path.join(WORK_DIR, file1))
        if file2:
            F2 = open(os.path.join(WORK_DIR, file2))

    line = F1.readline()
    sample_read_id = line.split(' ')[0]
    F1.seek(0)
    if sample_read_id.startswith('>'):
        F1.close()
        if file2:
            F2.close()
        read_data = studyFastaReads(read_data)
        return read_data
    
    read_ids_paired = True
    seqLen1 = 0
    seqLen2 = 0
    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    #minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    read_id_1 = None
    read_id_2 = None
    for i, line1 in enumerate(F1):
        if file2:
            line2 = F2.readline()
            if not line2:
                comment = "Number of reads differs between {} and {} after {}".format(file1, file2, readNumber)
                LOG.write(comment+"\n")
                read_data["problem"].append(comment)
                break
            if i % 4 == 0 and read_ids_paired:
                read_id_1 = line1.split(' ')[0] # get part up to first space, if any 
                read_id_2 = line2.split(' ')[0] # get part up to first space, if any 
                if not read_id_1 == read_id_2:
                    diff = findSingleDifference(read_id_1, read_id_2)
                    if diff == None or sorted((read_id_1[diff[0]:diff[1]], read_id_2[diff[0]:diff[1]])) != ('1', '2'):
                        read_ids_paired = False
                        read_data["problem"].append("id_mismatch at read %d: %s vs %s"%(readNumber+1, read_id_1, read_id_2))
            if i % 4 == 1:
                seqLen1 = len(line1)-1
                seqLen2 = len(line2)-1
                totalReadLength += (seqLen1 + seqLen2)
                maxReadLength = max(maxReadLength, seqLen1, seqLen2) 
                #minReadLength = min(minReadLength, seqLen1, seqLen2)
                readNumber += 1
            elif i % 4 == 3:
                minQualScore = min(minQualScore + line1.rstrip() + line2.rstrip())
                maxQualScore = max(maxQualScore + line1.rstrip() + line2.rstrip())
                if seqQualLenMatch:
                    if not (seqLen1 == len(line1)-1 and seqLen2 == len(line2)-1):
                        readId = [read_id_1, read_id_2][seqLen1 != len(line1)-1]
                        seqQualLenMatch = False
                        comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, readId)
                        if len(read_data["problem"]) < 50:
                            read_data["problem"].append(comment)
                        read_data["seq_qual_lengths_differ"] = True
                        LOG.write(comment+"\n")
        else: # no line2 -- single-end reads
            if i % 4 == 0:
                read_id = line1.split(' ')[0] # get part up to first space, if any 
            if i % 4 == 1:
                seqLen1 = len(line1)-1
                totalReadLength += seqLen1
                maxReadLength = max(maxReadLength, seqLen1) 
                # minReadLength = min(minReadLength, seqLen1)
                readNumber += 1
            elif i % 4 == 3:
                minQualScore = min(minQualScore + line1.rstrip())
                maxQualScore = max(maxQualScore + line1.rstrip())
                if seqQualLenMatch:
                    if not (seqLen1 == len(line1)-1):
                        seqQualLenMatch = False
                        comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, read_id)
                        if len(read_data["problem"]) < 50:
                            read_data["problem"].append(comment)
                        read_data["seq_qual_lengths_differ"] = True
                        LOG.write(comment+"\n")
        if False and readNumber % 100000 == 0:
            LOG.write("number of reads and bases tudied so far: \t{}\t{}\n".format(readNumber, totalReadLength))
            LOG.flush()

    F1.close()
    if file2:
        F2.close()

    avgReadLength = 0
    if readNumber:
        avgReadLength = totalReadLength/readNumber
    if file2:
        avgReadLength/=2
    read_data['avg_len'] = avgReadLength
    read_data['max_read_len'] = maxReadLength
    #read_data['min_read_len'] = minReadLength
    read_data['num_reads'] = readNumber
    read_data['num_bases'] = totalReadLength
    read_data['sample_read_id'] = sample_read_id 

    read_data['inferred_platform'] = inferPlatform(sample_read_id, maxReadLength)
    read_data['length_class'] = ["short", "long"][maxReadLength >= MAX_SHORT_READ_LENGTH]
    if file2 and maxReadLength >= MAX_SHORT_READ_LENGTH:
        comment = "paired reads appear to be long, expected short: %s"%read_data['delim'].join(read_data['files'])
        LOG.write(comment+"\n")
        read_data['problem'].append(comment)
    LOG.write("analysis of {} shows read_number = {} and total_bases = {}\n".format(":".join(read_data['files']), readNumber, totalReadLength))
    if False: #debugging
        comment = ""
        for key in read_data:
            comment += "study_reads: {} is {}\n".format(key, read_data[key])
        LOG.write(comment)
    LOG.write("duration of study_reads was %d seconds\n"%(time() - func_start))
    LOG.flush()
    return read_data

def studyFastaReads(read_data):
    """
    assume format is fasta
    count reads, calc total length, mean, max, min
    """
    func_start = time()
    seq = ""
    seqLen = 0
    totalReadLength = 0
    maxReadLength = 0
    #minReadLength = 1e6
    readNumber = 0
    file1 = read_data['files'][0]
    if file1.endswith("gz"):
        F = gzip.open(os.path.join(WORK_DIR, file1))
    else:
        F = open(os.path.join(WORK_DIR, file1))
    for line in F:
        if line.startswith(">"):
            readNumber += 1
            if seq:
                seqLen = len(seq)
                totalReadLength += seqLen
                maxReadLength = max(maxReadLength, seqLen) 
                #minReadLength = min(minReadLength, seqLen)
                seq = ""
            else:
                sample_read_id = line.split()[0]
        else:
            seq += line.rstrip()
    if seq:
        seqLen = len(seq)
        totalReadLength += seqLen
        maxReadLength = max(maxReadLength, seqLen) 
        #minReadLength = min(minReadLength, seqLen)

    avgReadLength = totalReadLength/readNumber
    read_data['avg_len'] = avgReadLength
    read_data['max_read_len'] = maxReadLength
    #read_data['min_read_len'] = minReadLength
    read_data['num_reads'] = readNumber
    read_data['sample_read_id'] = sample_read_id 
    read_data['platform'] = 'fasta'
    read_data['inferred_platform'] = inferPlatform(sample_read_id, maxReadLength)
    read_data['length_class'] = ["short", "long"][maxReadLength >= MAX_SHORT_READ_LENGTH]

    LOG.write("duration of studyFastaReads was %d seconds\n"%(time() - func_start))
    return read_data

def rewrite_reads(read_data, max_bases=0):
    """
    Read both files. Verify read ID are paired. Determine avg read length. Update details['reads'].
    """
    func_start = time()
    LOG.write("subsamplePairedReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start - START_TIME))

    output_read_data = read_data.copy()
    if 'num_bases' in read_data:
        output_read_data['num_bases_original'] = read_data['num_bases']
    if 'num_reads' in read_data:
        output_read_data['num_reads_original'] = read_data['num_reads']

    file1 = read_data['files'][0]
    file2 = None
    if len(read_data['files']) > 1:
        file2 = read_data['files'][1]
    suffix = "_subsampled.fastq"

    out_file1, out_file2 = file1, file2
    if file1.endswith("gz"):
        F1 = gzip.open(os.path.join(WORK_DIR, file1))
        out_file1 = file1[:-3]+suffix
        if file2:
            F2 = gzip.open(os.path.join(WORK_DIR, file2))
            out_file2 = file2[:-3]+suffix
    elif file1.endswith("bz2"):
        F1 = bz2.BZ2File(os.path.join(WORK_DIR, file1))
        out_file1 = file1[:-4]+suffix
        if file2:
            F2 = bz2.BZ2File(os.path.join(WORK_DIR, file2))
            out_file2 = file2[:-4]+suffix
    else:
        F1 = open(os.path.join(WORK_DIR, file1))
        out_file1 = file1+suffix
        if file2:
            F2 = open(os.path.join(WORK_DIR, file2))
            out_file2 = file2+suffix

    output_read_data['files']=[out_file1]
    if file2:
        output_read_data['files'].append(out_file2)

    output_read_data['files'] = (out_file1, out_file2)[:len(read_data['files'])]
    OF1 = open(os.path.join(WORK_DIR, out_file1), 'w')
    if file2:
        OF2 = open(os.path.join(WORK_DIR, out_file2), 'w')
    output_read_data['problem_reads']=[]
    reads_written = 0
    bases_written = 0
    record1 = ''
    record2 = ''
    problem = False
    for i, line1 in enumerate(F1):
        record1 += line1
        if file2:
            line2 = F2.readline()
            record2 += line2
        record_index = i % 4
        if record_index == 0:
            read1_id = line1.rstrip()
            if file2:
                read2_id = line2.rstrip()
        if record_index == 1:
            read1_length = len(line1) - 1
            if file2:
                read2_length = len(line2) - 1
        elif record_index == 3:
            qual1_length = len(line1) - 1
            if file2:
                qual2_length = len(line2) - 1
            if read1_length != qual1_length:
                problem = True
                if len(output_read_data['problem_reads']) < 50:
                    output_read_data['problem_reads'].append("Read {} from {}: read length = {}, qual length = {}".format(read1_id, file1, read1_length, qual1_length))
            if file2 and read2_length != qual2_length:
                problem = True
                if len(output_read_data['problem_reads']) < 50:
                    output_read_data['problem_reads'].append("Read {} from {}: read length = {}, qual length = {}".format(read2_id, file2, read2_length, qual2_length))

            if not problem:
                OF1.write(record1)
                bases_written += read1_length
                if file2:
                    OF2.write(record2)
                    bases_written += read2_length
                reads_written += 1
                if max_bases and bases_written >= max_bases:
                    break
            record1 = ''
            record2 = ''
            problem = False

    OF1.close()
    if file2:
        OF2.close()

    output_read_data['num_reads'] = reads_written
    output_read_data['num_bases'] = bases_written
    avgReadLength = bases_written/reads_written
    if file2:
        avgReadLength /= 2 # paired reads
    output_read_data['avg_len'] = avgReadLength

    LOG.write("duration of rewrite_reads was %d seconds\n"%(time() - func_start))
    return output_read_data

def inferPlatform(read_id, maxReadLength):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    if read_id.startswith(">"):
        return "fasta"
    if maxReadLength < MAX_SHORT_READ_LENGTH:
        # example illumina read id
        #@D00553R:173:HG53VBCXY:2:1101:1235:2074 1:N:0:ACAGTGAT
        parts = read_id.split(":")
        if len(parts) == 3:
            return "iontorrent"
        if len(parts) > 4:
            return "illumina"
        if re.match(r"@[A-Z]\S+:\d+:\S+:\d+:\d+:\d+:\d+ \S+:\S+:\S+:\S+$", read_id):
            return "illumina" # newer illumina
        if re.match(r"@\S+:\S+:\S+:\S+:\S+#\S+/\S+$", read_id):
            return "illumina" # older illumina
        if re.match(r"@[^:]+:[^:]+:[^:]+$", read_id):
            return "iontorrent" # 
        if re.match(r"@[SED]RR\d+\.\d+", read_id):
            return "illumina" # default short fastq type 
    else: # one of the long read types (pacbio or nanopore)
    # todo: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
        if re.match(r"@\S+/\S+/\S+_\S+$", read_id): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
            return "pacbio" # 
        if re.match(r"@\S+/\S+$", read_id): #@<MovieName>/<ZMW_number> 
            return "pacbio" # 
    #@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
        if re.match(r"@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", read_id): #based on one example, need to test more 
            return "nanopore" # 
        return "pacbio" # default long fastq type
    # if we get here, we failed to recognize what read type it is, default to illumina
    LOG.write("inferPlatform defaulting to 'illumina'\n")
    return "illumina"

def trimGalore(details, threads=1):
    startTrimTime = time()
    LOG.write("\ntrimGalore() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    command = ["trim_galore", "--version"]
    proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE)
    version_text = proc.stdout.read().decode().strip()
    proc.wait()
    m = re.search(r"(version\s+\S+)", version_text)
    if m:
        trim_galore_version = "trim_galore " + m.group(1)
        details["version"]['trim_galore'] = trim_galore_version
    if "trim report" not in details:
        details["trim report"] = {}
    toRegister = {} # save trimmed reads to register after iterating dictionary to avoid error "dictionary changed size during iteration"
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'short' and not details['reads'][reads]['platform'] == 'fasta':
            command = ['trim_galore', '-j', str(min(threads, 8)), '-o', '.']
            if len(details['reads'][reads]['files']) > 1:
                command.extend(["--paired", details['reads'][reads]['files'][0], details['reads'][reads]['files'][1]])
            else:
                command.append(reads)

            LOG.write("command: "+" ".join(command)+"\n")
            proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
            trimGaloreStderr = proc.stderr.read()
            return_code = proc.wait()
            LOG.write("return code = %d\n"%return_code)
            trimReads = re.findall(r"Writing validated paired-end read \d reads to (\S+)", trimGaloreStderr)
            if not trimReads:
                trimReads = re.findall(r"Writing final adapter and quality trimmed output to (\S+)", trimGaloreStderr)
            LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
            if not trimReads:
                comment = "trim_galore did not name trimmed reads output files in stderr"
                LOG.write(comment+"\n")
                details['reads'][reads]['problem'].append(comment)
                continue
            comment = "trim_galore, input %s, output %s"%(reads, ":".join(trimReads))
            LOG.write(comment+"\n")
            trimmed_reads = trimReads[0]
            if len(trimReads) > 1:
                trimmed_reads += ":"+trimReads[1]
            details['reads'][reads]['trimmed_reads'] = trimmed_reads
            details["pre-assembly transformation"].append(comment)
            toRegister[trimmed_reads] = reads

            trimReports = re.findall(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
            LOG.write("re.findall for trim reports returned %s\n"%str(trimReports))
            details["trim report"][reads]=[]
            for f in trimReports:
                filename_base = os.path.basename(f)
                shutil.move(f, os.path.join(DETAILS_DIR, filename_base))
                details["trim report"][reads].append(filename_base)

            
    for trimReads in toRegister:
        registerReads(trimReads, details, supercedes=toRegister[trimReads], max_bases=details['max_bases'])

    LOG.write("trim_galore trimming completed, duration = %d seconds\n\n\n"%(time()-startTrimTime))

def sampleReads(filename, details=None):
    srf_time = time()
    LOG.write("sampleReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(srf_time)), srf_time-START_TIME))
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
    text = F.read(Default_bytes_to_sample) #read X number of bytes for text sample
    F.close()

    LOG.write("  file text sample %s:\n%s\n\n"%(filename, text[0:50]))
    lines = text.split("\n")
    readLengths = []
    if len(lines) < 2:
        comment = "in sampleReads for %s: text sample (length %d) lacks at least 2 lines"%(filename, len(text))
        LOG.write(comment+"\n")
        details["problem"].append(comment)
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
        LOG.write(comment+"\n")
        details["problem"].append(comment)
    if len(read_id_sample) > 1:
        read_id_sample = read_id_sample[:-1] # last entry might be truncated, avoid it
    LOG.write("read type %s, maximum read length %.1f\n"%(read_format, max_read_length))
    return read_id_sample, max_read_length

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

def categorize_anonymous_read_files(args, details):
    LOG.write("categorize_anonymous_read_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("  files=%s\n"%("\t".join(args.anonymous_reads)))

    nonSraFiles = []
    # first pull out any SRA files for special treatment
    for item in args.anonymous_reads:
        m = re.match(r"^([SED]RR\d+)$", item)
        if m:
            sra = m.group(1)
            if sra not in args.sra:
                args.sra.append(sra)
        else:
            nonSraFiles.append(item)

    # now proceed with any non-sra files
    read_file_type = {}
    read_id_sample = {}
    singleFiles = []
    pairedFiles = []
    for item in nonSraFiles:
        if ":" in item or "%" in item:
            pairedFiles.append(item)
        else:
            read_id_sample[item], max_read_length = sampleReads(item, details)
            read_file_type[item] = inferPlatform(read_id_sample[item][0], max_read_length)
            comment = "interpreting %s type as %s"%(item, read_file_type[item])
            LOG.write(comment+"\n")
            #details["pre-assembly transformation"].append(comment)
            if read_file_type[item] is not None:
                singleFiles.append(item)

    # try to find paired files
    membersOfPairs = set()
    for i, filename1 in enumerate(singleFiles[:-1]):
        for filename2 in singleFiles[i+1:]:
            singleDiff = findSingleDifference(filename1, filename2)
# singleDiff will be not None if the strings match at all but one character (presumably '1' vs '2')
            if singleDiff and singleDiff[0] > 0 and singleDiff[1]-singleDiff[0] == 1:
                charBefore = filename1[singleDiff[0]-1]
                if charBefore.isdigit():
                    continue # changes in multi-digit numbers are not indicative of paired reads
                diffChars = (filename1[singleDiff[0]], filename2[singleDiff[0]])
                pair = None
                if diffChars[0] == '1' and diffChars[1] == '2':
                    pair = (filename1, filename2)
                elif diffChars[1] == '1' and diffChars[0] == '2':
                    pair = (filename2, filename1)
                if pair:
                    comment = "candidate paired files: %s  %s"%pair
                    LOG.write(comment+"\n")
                    details["problem"].append(comment)
                    if read_file_type[filename1] != read_file_type[filename2]:
                        comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
                        LOG.write(comment+"\n")
                        details["problem"].append(comment)
                        continue
                    pairedFiles.append(pair[0] + ":" + pair[1])

    # now go over all pairs to test for matching types and matching read IDs
    valid_pairs = set()
    valid_singles = set()
    for item in pairedFiles:
        if ":" in item:
            filename1, filename2 = item.split(":") 
        elif "%" in item:
            filename1, filename2 = item.split("%") 
        else:
            comment = "failed to find ':' or '%%' in file pair: %s" % (item)
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            continue
        if filename1 not in read_file_type:
            read_id_sample[filename1], max_read_length = sampleReads(filename1, details)
            read_file_type[filename1] = inferPlatform(read_id_sample[filename1][0], max_read_length)
            comment = "interpreting %s type as %s"%(filename1, read_file_type[filename1])
            LOG.write(comment+"\n")
            #details["pre-assembly transformation"].append(comment)
        if filename2 not in read_file_type:
            read_id_sample[filename2], max_read_length = sampleReads(filename2, details)
            read_file_type[filename2] = inferPlatform(read_id_sample[filename2][0], max_read_length)
            comment = "interpreting %s type as %s"%(filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            #details["pre-assembly transformation"].append(comment)

        read_types_match = True
        # test if read types are the same
        if read_file_type[filename1] != read_file_type[filename2]:
            comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            #details["problem"].append(comment)
            read_types_match = False
        read_file_type[item] = read_file_type[filename1] #easier to retrieve later

        # test if read IDs match between files
        ids_paired = True
        for idpair in zip(read_id_sample[filename1], read_id_sample[filename2]):
            if idpair[0] == idpair[1]:
                continue
            if idpair[0][:-1] == idpair[1][:-1]:
                continue
            ids_paired = False
            comment = "Read IDs do not match for %s(%s) vs %s(%s)"%(filename1, idpair[0], filename2, idpair[1])
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            singleFiles.extend((filename1, filename2)) #move over to single files
            break
        if read_types_match and ids_paired:
            valid_pairs.add(item)
            membersOfPairs.add(filename1)
            membersOfPairs.add(filename2)
        else: #move over to single files
            valid_singles.add(filename1)
            valid_singles.add(filename2)
    
    # some items on singleFiles may not be valid (may be paired up)
    for item in singleFiles:
        if item not in membersOfPairs:
            valid_singles.add(item)

    for item in valid_pairs.union(valid_singles):
        if read_file_type[item] == "illumina":
            args.illumina.append(item)
        elif read_file_type[item] == "iontorrent":
            args.iontorrent.append(item)
        elif read_file_type[item] == "pacbio":
            args.pacbio.append(item)
        elif read_file_type[item] == "nanopore":
            args.nanopore.append(item)
        elif read_file_type[item] == "fasta":
            if not hasattr(args, 'fasta'):
                setattr(args, 'fasta', [])
            args.fasta.append(item)
        else:
            comment = "Cannot decide read type for item "+item
            LOG.write(comment+"\n")
            details['problem'].append(comment)
        registerReads(item, details, platform=read_file_type[item], interleaved = (args.interleaved and item in args.interleaved), max_bases=args.max_bases)

    return

def get_sra_runinfo(run_accession, log=None):
    """ take sra run accession (like SRR123456)
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    Altered from version in sra_tools to handle case of multiple sra runs returned by query.
    If efetch doesn't work, try scraping the web page.
    """
    LOG.write("get_sra_runinfo(%s)\n"%run_accession)
    if run_accession.endswith(".sra"):
        run_accession = run_accession[:-4]
    runinfo = None
    runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="+run_accession
    text = urllib2.urlopen(runinfo_url).read()
    if text.startswith("Run"):
        lines = text.split("\n")
        keys   = lines[0].split(",")
        for line in lines[1:]:  # there might be multiple rows, only one of which is for our sra run accession
            if line.startswith(run_accession):
                values = line.split(",")
                runinfo = dict(zip(keys, values))
    if runinfo:
        if runinfo['Platform'] not in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
            # rescue case of mis-alignment between keys and values
            if log:
                log.write("problem in get_sra_runinfo: sra.cgi returned:\n"+text+"\n")
            for val in values:
                if val in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
                    runinfo['Platform'] = val
                    break
        if not runinfo['LibraryLayout'].startswith(('PAIRED', 'SINGLE')):        
            if log:
                log.write("Need to search for LibraryLayout: bad value: %s\n"%runinfo['LibraryLayout'])
            for val in values:
                if val.startswith(('PAIRED', 'SINGLE')):
                    runinfo['LibraryLayout'] = val
                    break

    if not runinfo:
        if log:
            log.write("Problem, normal runinfo request failed. Trying alternative from web page.\n")
        # screen-scrape
        runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run="+run_accession
        text = urllib2.urlopen(runinfo_url).read()
        runinfo = {}
        if re.search("<td>Illumina", text, re.IGNORECASE):
            runinfo['Platform'] = 'ILLUMINA'
        elif re.search("<td>PacBio", text, re.IGNORECASE):
            runinfo['Platform'] = 'PACBIO_SMRT'
        elif re.search("<td>Oxford", text, re.IGNORECASE):
            runinfo['Platform'] = 'OXFORD_NANOPORE'
        elif re.search("<td>Ion Torrent", text, re.IGNORECASE):
            runinfo['Platform'] = 'ION_TORRENT'

        if re.search("<td>SINGLE</td>", text, re.IGNORECASE):
            runinfo['LibraryLayout'] = 'SINGLE'
        elif re.search("<td>PAIRED", text, re.IGNORECASE):
            runinfo['LibraryLayout'] = 'PAIRED'
    return runinfo

def fetch_one_sra(sra, run_info=None, log=sys.stderr, usePrefetch=False):
    """ requires run_info to know which program to use
    """
    if not run_info:
        run_info = get_sra_runinfo(sra, log)

    if usePrefetch:
        command = ["prefetch", sra]
        log.write("command = "+" ".join(command)+"\n")
        return_code = subprocess.call(command, shell=False, stderr=log)
        log.write("return_code = %d\n"%(return_code))

    command = ["fasterq-dump", "--split-files", sra] # but not appropriate for pacbio or nanopore
    if run_info['Platform'].startswith("PACBIO") or run_info['Platform'].startswith("OXFORD_NANOPORE"):
        command = ['fastq-dump', sra]
    stime = time()
    log.write("command = "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=log)
    log.write("return_code = %d\n"%(return_code))
    if return_code != 0:
        log.write("Problem, return code was %d\n"%(return_code))

        log.write("Try one more time.\n")
        return_code = subprocess.call(command, shell=False, stderr=LOG)
        log.write("Return code on second try was %d\n"%return_code)
        if return_code != 0:
            LOG.write("Giving up on %s\n"%sra)
    LOG.write("fetch_one_sra time=%d seconds\n"%(time()-stime))
    return

def fetch_sra_files(sra_ids, details):
    """ 
    fetch each sra item and register the reads
    """
    LOG.write("fetch_sra_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("sra_ids="+" ".join(sra_ids)+"\n")
    for sra in sra_ids:
        sraFull = sra
        sra = sraFull.replace(".sra", "")
        runinfo = get_sra_runinfo(sra)
        if not runinfo:
            LOG.write("runinfo for %s was empty, giving up\n"%sra)
            continue
        LOG.write("Runinfo for %s reports platform = %s and LibraryLayout = %s\n"%(sra, runinfo["Platform"], runinfo['LibraryLayout']))

        fetch_one_sra(sraFull, run_info=runinfo, log=LOG)
        fastqFiles = glob.glob(sra+"*fastq")
        LOG.write("Fastq files from sra: %s\n"%str(fastqFiles))
        processSraFastqFiles(fastqFiles, details, runinfo)
    return

def processSraFastqFiles(fastqFiles, details, run_info=None):
    """ manipulate multiple (or single) fastq files from one SRA runId and register reads """
    comment = "processSraFastqFiles(%s)"%",".join(fastqFiles)
    details['problem'].append(comment)
    LOG.write(comment+"\n")
    item = None
    m = re.match(r"([SED]RR\d+)", os.path.basename(fastqFiles[0]))
    if not m:
        comment = "supposed sra fastq file does not start with [SED]RRnnnn"
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return
    sra = m.group(1)
    for fq in fastqFiles:
        if not os.path.basename(fq).startswith(sra):
            comment = "Problem: not all fastqFiles passed to processSraFastqFiles() begin with %s: %s"%(sra, ",".join(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            return
    if not run_info:
        run_info = get_sra_runinfo(sra)

    if run_info['LibraryLayout'].startswith("PAIRED"):
        if len(fastqFiles) == 2:
            item = ":".join(sorted(fastqFiles)[:2])
            comment = "runinfo[LibraryLayout] == PAIRED: item = %s"%item
            details['problem'].append(comment)
            LOG.write(comment+"\n")
        else:
            comment = "for PAIRED library %s, number of files was %s, expected 2: %s"%(sra, len(fastqFiles), str(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            if len(fastqFiles) == 1:
                item = fastqFiles[0] # interpret as single-end, perhaps an SRA metadata mistake
                comment = "interpret library %s as single-end"%sra
                details['problem'].append(comment)
                LOG.write(comment+"\n")
    if not item: # library layout single or failed above
        if len(fastqFiles) == 1:
            item = fastqFiles[0]
            comment = "runinfo[LibraryLayout] == %s: item = %s"%(run_info['LibraryLayout'], item)
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            
        elif len(fastqFiles) > 1:
            comment = "LibraryLayout=%s; Platform=%s: multiple files = %s"%(run_info['LibraryLayout'], run_info['Platform'], ",".join(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")

            concatenate_command = "cat %s > %s/%s.fastq"%(" ".join(fastqFiles), WORK_DIR, sra)
            LOG.write("concatenate command:"+concatenate_command+"\n")
            subprocess.call(concatenate_command, shell=True)
            item = sra+".fastq"
            comment = "for library %s, list of files was %s, concatenated to %s"%(sra, str(fastqFiles), item)
            details['problem'].append(comment)
            LOG.write(comment+"\n")
        if not item:
            comment = "for %s no fastq file found"%sra
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            # failed on that sra
   
    platform = None
    if run_info["Platform"] == "ILLUMINA":
        platform = "illumina"
    elif run_info["Platform"] == "ION_TORRENT":
        platform = "iontorrent"
    elif run_info["Platform"] == "PACBIO_SMRT":
        platform = "pacbio"
    elif run_info["Platform"] == "OXFORD_NANOPORE":
        platform = "nanopore"
    if not platform:
        ids, max_length = sampleReads(fastqFiles[0], details)
        platform = ["illumina", "pacbio"][max_length > MAX_SHORT_READ_LENGTH]
    registerReads(item, details, platform=platform, max_bases=details['max_bases'])
    return

def writeSpadesYamlFile(details):
    LOG.write("writeSpadesYamlFile: elapsed seconds = %f\n"%(time()-START_TIME))
    outfileName = "spades_yaml_file.txt"
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    for platform in details['platform']:
        LOG.write(platform+": "+", ".join(details['platform'][platform])+"\n")
    
    single_end_reads = []
    paired_end_reads = [[], []]
    mate_pair_reads = [[], []]
    interleaved_reads = []

    shortReadItems = []
    if 'illumina' in details['platform']:
        shortReadItems.extend(details['platform']['illumina'])
    if 'iontorrent' in details['platform']:
        shortReadItems.extend(details['platform']['iontorrent'])
    if 'fasta' in details['platform']:
        shortReadItems.extend(details['platform']['fasta'])
    for item in shortReadItems:
        LOG.write("process item {}\n".format(item))
        if ":" in item:
            f = details['reads'][item]['files'][0]
            paired_end_reads[0].append(f)
            f = details['reads'][item]['files'][1]
            paired_end_reads[1].append(f)
        elif "%" in item:
            f = details['reads'][item]['files'][0]
            mate_pair_reads[0].append(f)
            f = details['reads'][item]['files'][1]
            mate_pair_reads[1].append(f)
        else:
            f = details['reads'][item]['files'][0]
            if 'interleaved' in details['reads'][item]:
                interleaved_reads.append(f)
            else:
                single_end_reads.append(f)

    precedingElement=False
    if single_end_reads:
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(single_end_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if interleaved_reads:
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"paired-end\",\n    interlaced reads: [\n        \"")
        OUT.write("\",\n        \"".join(interleaved_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if paired_end_reads[0]:
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"fr\",\n")
        OUT.write("    type: \"paired-end\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(paired_end_reads[0]))
        OUT.write("\"\n    ],\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(paired_end_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement = True
    if mate_pair_reads[0]:
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"rf\",\n")
        OUT.write("    type: \"mate-pairs\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[0]))
        OUT.write("\"\n    ]\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement = True
    if details['platform']['pacbio']:
        pacbio_reads = []
        for item in details['platform']['pacbio']:
            f = details['reads'][item]['files'][0]
            pacbio_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if details['platform']['nanopore']:
        nanopore_reads = []
        for item in details['platform']['nanopore']:
            f = details['reads'][item]['files'][0]
            nanopore_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if False and details['platform']['fasta']:
        fasta_reads = []
        for item in details['platform']['fasta']:
            f = details['reads'][item]['files'][0]
            fasta_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(fasta_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True

    OUT.write("]\n")
    OUT.close()
    return(outfileName)    

def runQuast(contigsFile, args, details):
    LOG.write("runQuast() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    quastDir = "quast_out"
    quastCommand = ["quast.py",
                    "-o", quastDir,
                    "-t", str(args.threads),
                    "--min-contig", str(args.min_contig_length),
                    contigsFile]
    LOG.write("running quast: "+" ".join(quastCommand)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(quastCommand, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("return code = %d\n"%return_code)
    if return_code == 0:
        shutil.move(os.path.join(quastDir, "report.html"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.html"))
        shutil.move(os.path.join(quastDir, "report.tsv"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.tsv"))
        shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.txt"))
        shutil.move(os.path.join(quastDir, "transposed_report.txt"), os.path.join(DETAILS_DIR, args.prefix+"quast_transposed_report.txt"))
        shutil.move(os.path.join(quastDir, "transposed_report.tsv"), os.path.join(DETAILS_DIR, args.prefix+"quast_transposed_report.tsv"))
        details["quast_transposed_txt"] = "details/"+args.prefix+"quast_transposed_report.txt"
        details["quast_transposed_tsv"] = "details/"+args.prefix+"quast_transposed_report.tsv"
        details["quast_txt"] = "details/"+args.prefix+"quast_report.txt"
        details["quast_tsv"] = "details/"+args.prefix+"quast_report.tsv"
        details["quast_html"] = "details/"+args.prefix+"quast_report.html"
        quastCommand = ["quast.py", "--version"]
        proc = subprocess.Popen(quastCommand, shell=False, stdout=subprocess.PIPE)
        version_text = proc.stdout.read().decode()
        proc.wait()
        details["version"]["quast"] = version_text

def filterContigsByMinLength(inputContigs, details, min_contig_length=300, min_contig_coverage=5, threads=1, prefix=""):
    """ 
    Write only sequences at or above min_length to output file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("filterContigsByMinLength(%s) \n"%(inputContigs))
    report = {}
    shortReadDepth = None
    longReadDepth = None
    bamFiles = []
    totalShortReadBases = 0
    totalLongReadBases = 0
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'short':
            bam = runBowtie(inputContigs, reads, details, threads=threads, outformat='bam')
            if bam:
                bamFiles.append(bam)
            if 'num_bases' in details['reads'][reads]:
                totalShortReadBases += details['reads'][reads]['num_bases']
    if bamFiles:
        shortReadDepth = calcReadDepth(bamFiles, details)
        if shortReadDepth and 'average_coverage' in shortReadDepth:
            report['average depth (short reads)'] = shortReadDepth['average_coverage']
    bamFiles = []
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'long':
            bam = runMinimap(inputContigs, reads, details, threads=threads, outformat='bam')
            if bam:
                bamFiles.append(bam)
            if 'num_bases' in details['reads'][reads]:
                totalLongReadBases += details['reads'][reads]['num_bases']
    if bamFiles:
        longReadDepth = calcReadDepth(bamFiles, details)
        if longReadDepth and 'average_coverage' in longReadDepth:
            report['average depth (long reads)'] = longReadDepth['average_coverage']
    report['min_contig_length_threshold'] = "%d"%min_contig_length
    report["min_contig_coverage_threshold"] = "%.1f"%min_contig_coverage
    if totalShortReadBases:
        report['total_short_read_bases'] = totalShortReadBases
        LOG.write("Total short read bases = %d\t"%totalShortReadBases)
    if totalLongReadBases:
        report['total_long_read_bases'] = totalLongReadBases
        LOG.write("Total long read bases = %d\t"%totalLongReadBases)
    num_good_contigs = num_bad_contigs = 0
    suboptimalContigs = ""
    total_seq_length = 0
    weighted_short_read_coverage = 0
    weighted_long_read_coverage = 0
    outputContigs = re.sub(r"\..*", "_depth_cov_filtered.fasta", inputContigs)
    LOG.write("writing filtered contigs to %s\n"%outputContigs)
    with open(inputContigs) as IN:
        with open(outputContigs, 'w') as OUT:
            seqId=None
            seq = ""
            contigIndex = 1
            num_circular_contigs = 0
            line = "1"
            while line:
                line = IN.readline()
                m = re.match(r">(\S+)", line)
                if m or not line: 
                    if seq:
                        contigId = ">"+prefix+"contig_%d"%contigIndex
                        contigInfo = " length %5d"%len(seq)
                        if len(seq) < min_contig_length:
                            suboptimalContigs += contigId+contigInfo+"\n"+seq+"\n"
                            num_bad_contigs += 1
                            continue
                        contigIndex += 1
                        short_read_coverage = 0
                        long_read_coverage = 0
                        passes_coverage_threshold = not (report['short read coverage available'] and report['long read coverage available'])
                        if shortReadDepth and seqId in shortReadDepth:
                            short_read_coverage, normalizedDepth = shortReadDepth[seqId]
                            contigInfo += " coverage %.01f normalized_cov %.2f"%(short_read_coverage, normalizedDepth)
                            passes_coverage_threshold = short_read_coverage >= min_contig_coverage
                        if longReadDepth and seqId in longReadDepth:
                            long_read_coverage, normalizedDepth = longReadDepth[seqId]
                            contigInfo += " longread_coverage %.01f normalized_longread_cov %.2f"%(long_read_coverage, normalizedDepth)
                            passes_coverage_threshold = long_read_coverage >= min_contig_coverage
                        if passes_coverage_threshold:
                            OUT.write(contigId+contigInfo+"\n")
                            for i in range(0, len(seq), 60):
                                OUT.write(seq[i:i+60]+"\n")
                            num_good_contigs += 1
                            if short_read_coverage:
                                weighted_short_read_coverage += short_read_coverage * len(seq)
                            if long_read_coverage:
                                weighted_long_read_coverage += long_read_coverage * len(seq)
                            total_seq_length += len(seq)
                        else:
                            suboptimalContigs += contigId+contigInfo+"\n"+seq+"\n"
                            num_bad_contigs += 1
                        if "circular=true" in line:
                            num_circular_contigs += 1
                        seq = ""
                    if m:
                        seqId = m.group(1)
                elif line:
                    seq += line.rstrip()
    if total_seq_length:
        if weighted_short_read_coverage:
            report['average short read coverage'] = "%.3f"%(weighted_short_read_coverage / total_seq_length)
        if weighted_long_read_coverage:
            report['average long read coverage'] = "%.3f"%(weighted_long_read_coverage / total_seq_length)
    report['num contigs above thresholds'] = "%d"%num_good_contigs
    report['num contigs below thresholds'] = "%d"%num_bad_contigs
    report['total length of good contigs'] = "%d"%total_seq_length
    if num_circular_contigs:
        report['contigs predicted circular'] = "%d"%num_circular_contigs
    if suboptimalContigs:
        suboptimalContigsFile = "contigs_below_length_coverage_threshold.fasta"
        report["suboptimal contigs file"] = suboptimalContigsFile
        suboptimalContigsFile = os.path.join(DETAILS_DIR, suboptimalContigsFile)
        with open(suboptimalContigsFile, "w") as SUBOPT:
            SUBOPT.write(suboptimalContigs)
    if os.path.getsize(outputContigs) < 10:
        LOG.write("failed to generate outputContigs, return None\n")
        return None
    details['contig_filtering'] = report
    comment = "filterContigsByMinLength, input %s, output %s"%(inputContigs, outputContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return outputContigs

def runBandage(gfaFile, details):
    imageFormat = ".svg"
    retval = None
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write("Bandage command =\n"+" ".join(command)+"\n")
        try:
            with open(os.devnull, 'w') as FNULL:
                return_code = subprocess.call(command, shell=False, stderr=FNULL)
            LOG.write("return code = %d\n"%return_code)
            if return_code == 0:
                retval = plotFile
                proc = subprocess.Popen(["Bandage", "--version"], shell=False, stdout=subprocess.PIPE)
                version_text = proc.stdout.read().decode().strip()
                proc.wait()
                details["Bandage"] = {
                        "version" : "Bandage "+version_text,
                    "plot" : plotFile,
                    "command line": command
                    }
                details["version"]["Bandage"] = version_text
            else:
                LOG.write("Error creating Bandage plot\n")
        except OSError as ose:
            comment = "Problem running Bandage: "+str(ose)
            LOG.write(comment+"\n")
            details['problem'].append(comment)
    return retval

def runUnicycler(details, threads=1, min_contig_length=0, prefix="", spades_exec=None):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("runUnicycler\n")
    proc = subprocess.Popen(["unicycler", "--version"], shell=False, stdout=subprocess.PIPE)
    version_text = proc.stdout.read().decode()
    if version_text:
        m = re.search(r"Unicycler\s+\S+", version_text, flags=re.IGNORECASE)
        if m:
            version_text = m.group(0) # entire match
    proc.wait()
    details["assembly"]['assembler'] = 'unicycler'
    details["assembly"]['version'] = version_text
    details["version"]["unicycler"] = version_text
    command = ["unicycler", "-t", str(threads), "-o", '.']
    if min_contig_length:
        command.extend(("--min_fasta_length", str(min_contig_length)))
    command.extend(("--keep", "2")) # keep files needed for re-run if necessary
    command.append("--no_pilon")  # we will run our own, if requested
    if spades_exec:
        command.extend("--spades_path", spades_exec);

    # put all read files on command line, let Unicycler figure out which type each is
    # apparently unicycler can only accept one read set in each class (I tried multiple ways to submit 2 paired-end sets, failed)
    short1 = None
    short2 = None
    unpaired = None
    long_reads = None
    for item in details['reads']:
        if 'superceded_by' in details['reads'][item]:
            continue
        files = details['reads'][item]['files'] 
        if details['reads'][item]['length_class'] == 'short':
            if len(files) > 1:
                short1 = files[0]
                short2 = files[1]
            else:
                unpaired = files[0]
        else:
            long_reads = files[0]

    if short1:
        command.extend(("--short1", short1, "--short2", short2))
    if unpaired:
        command.extend(("--unpaired", unpaired))
    if long_reads:
        command.extend(("--long", long_reads))
    # it is not quite right to send iontorrent data to spades through unicycler because the --iontorrent flag to spades will not be set

    LOG.write("Unicycler command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    LOG.flush()
    unicyclerStartTime = time()
    details["assembly"]['command_line'] = " ".join(command)
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big and unicycle.log is better
        return_code = subprocess.call(command, shell=False, stdout=FNULL)
    LOG.write("return code = %d\n"%return_code)

    if not (os.path.exists("assembly.fasta") and os.path.getsize("assembly.fasta")):
        comment = "First run of Unicycler resulted in no assembly, try again with more lenient parameters."
        LOG.write(comment+"\n")
        if 'problem' not in details['assembly']:
            details["assembly"]["problem"] = []
        details["assembly"]["problem"].append(comment)
        command.extend(("--mode", "bold", "--min_component_size", "300", "--min_dead_end_size", "300", "--depth_filter", "0.1"))
        comment = "re-run unicycler with command = "+" ".join(command)
        LOG.write(comment+"\n")
        details["assembly"]["problem"].append(comment)
        with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big and unicycle.log is better
            return_code = subprocess.call(command, shell=False, stdout=FNULL)
        LOG.write("return code = %d\n"%return_code)

    unicyclerEndTime = time()
    elapsedTime = unicyclerEndTime - unicyclerStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"]['assembly_time'] = elapsedHumanReadable
    details["assembly"]['assembly_seconds'] = elapsedTime

    LOG.write("Duration of Unicycler run was %s\n"%(elapsedHumanReadable))

    unicyclerLogFile = "unicycler.log"
    if os.path.exists("unicycler.log"):
        unicyclerLogFile = prefix+"unicycler.log"
        shutil.move("unicycler.log", os.path.join(DETAILS_DIR, unicyclerLogFile))

    if not os.path.exists("assembly.fasta"):
        comment = "Unicycler failed to generate assembly file. Check "+unicyclerLogFile
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        if 'problem' not in details['assembly']:
            details['assembly']['problem'] = []
        details["assembly"]["problem"].append(comment)
        return None

    assemblyGraphFile = prefix+"assembly_graph.gfa"
    shutil.move("assembly.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))

    contigsFile = "contigs.fasta"
    shutil.move("assembly.fasta", contigsFile) #rename to canonical name
    details["assembly"]["contigs.fasta file size"] = os.path.getsize(contigsFile)
    return contigsFile

def runSpades(details, args):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("runSpades: elapsed seconds = %f\n"%(time()-START_TIME))
    if args.illumina and args.iontorrent:
        comment = "SPAdes is not meant to process both Illumina and IonTorrent reads in the same run"
        if 'problem' not in details['assembly']:
            details['assembly']['problem'] = []
        details["assembly"]["problem"].append(comment)
        LOG.write(comment+"\n")

    command = ["spades.py", "--version"]
    proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    version_text = proc.stdout.read().decode().strip()
    if not version_text: # as of v3.14.1 this output came to stderr instead of stdout
        version_text = proc.stderr.read().decode().strip()
    proc.wait()
    details["assembly"]['assembler'] = 'SPAdes'
    details["assembly"]['version'] = version_text
    details["version"]["spades.py"] = version_text


    command = ["spades.py", "--threads", str(args.threads), "-o", "."]
    if args.recipe == 'single-cell':
        command.append("--sc")
    if args.iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(details)
    command.extend(["--dataset", yamlFile])
    if args.trusted_contigs:
        command.extend(["--trusted-contigs", args.trusted_contigs])
    if args.untrusted_contigs:
        command.extend(["--untrusted-contigs", args.untrusted_contigs])
    if hasattr(args, 'fasta') and args.fasta:
        command.append("--only-assembler")
    if args.memory:
        command.extend(["-m", str(args.memory)])
    if args.recipe == "meta-spades":
        command.append("--meta")
    if args.recipe == "rna-spades":
        command.append("--rna")
        #
        # Validate arguments for metagenomic spades. It can only run with
        # a single paired-end library.
        #
        #if len(single_end_reads) > 0 or len(paired_end_reads[0]) > 1:
        #    sys.stderr.write("SPAdes in metagenomics mode can only process a single paired-end read file\n")
        #    sys.exit(1);
    if args.recipe == "plasmid-spades":
        command.append("--plasmid")
    LOG.write("SPAdes command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    LOG.flush()
    spadesStartTime = time()

    details['assembly']['command_line'] = " ".join(command)
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("return code = %d\n"%return_code)

    contigsFile = "contigs.fasta"
    if return_code and not os.path.exists(contigsFile):
        comment = "spades return code = %d, attempt re-running without error correction"%return_code
        LOG.write(comment+"\n")
        if 'problem' not in details['assembly']:
            details['assembly']['problem'] = []
        details['assembly']['problem'].append(comment)

        #first run trimming unless it has already been run
        if not args.trim:
            trimGalore(details, threads=args.threads)
        yamlFile = writeSpadesYamlFile(details)
        command.append("--only-assembler")

        comment = "rerun SPAdes: command = "+" ".join(command)+"\n"
        LOG.write(comment+"\n")
        #details['assembly']['problem'].append(comment)
        details['assembly']['command_line'] = " ".join(command)

        with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
            return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
        LOG.write("return code = %d\n"%return_code)

    spadesEndTime = time()
    elapsedTime = spadesEndTime - spadesStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"]['assembly_time'] = elapsedHumanReadable

    LOG.write("Duration of SPAdes run was %s\n"%(elapsedHumanReadable))
    spadesLogFile = args.prefix+"spades.log"
    try:
        shutil.move("spades.log", os.path.join(DETAILS_DIR, spadesLogFile))
        assemblyGraphFile = args.prefix+"assembly_graph.gfa"
        shutil.move("assembly_graph_with_scaffolds.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))
    except Exception as e:
        LOG.write(str(e))
    if not os.path.exists(contigsFile):
        comment = "SPAdes failed to generate contigs file. Check "+spadesLogFile
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        details["problem"].append(comment)
        return None
    details["assembly"]["contigs.fasta size:"] = os.path.getsize(contigsFile)
    return contigsFile

def runMinimap(contigFile, longReadFastq, details, threads=1, outformat='sam'):
    #LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    """
    Map long reads to contigs by minimap2 (read paf-file, readsFile; generate paf file).
    """
    LOG.write("runMinimap(%s, %s, details, %d, %s)\n"%(contigFile, longReadFastq, threads, outformat))
    if 'minimap2' not in details['version']:
        command = ["minimap2", "--version"]
        proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE)
        proc.wait()
        version_text = proc.stdout.read().decode()
        details["version"]['minimap2'] = version_text.strip()
    # index contig sequences
    contigIndex = contigFile.replace(".fasta", ".mmi")
    command = ["minimap2", "-t", str(threads), "-d", contigIndex, contigFile] 
    tempTime = time() 
    LOG.write("minimap2 index command:\n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("minimap2 index return code = %d, time = %d seconds\n"%(return_code, time() - tempTime))
    if return_code != 0:
        return None

    # map long reads to contigs
    contigSam = contigFile.replace(".fasta", ".sam")
    command = ["minimap2", "-t", str(threads), "-a", "-o", contigSam, contigFile, longReadFastq]
    tempTime = time()
    LOG.write("minimap2 map command:\n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL)
    LOG.write("minimap2 map return code = %d, time = %d seconds\n"%(return_code, time() - tempTime))
    if return_code != 0:
        return None

    if outformat == 'sam':
        LOG.write('runMinimap returning %s\n'%contigSam)
        return contigSam

    else:
        contigsBam = convertSamToBam(contigSam, details, threads=threads)
        LOG.write('runMinimap returning %s\n'%contigsBam)
        return contigsBam
            

def convertSamToBam(samFile, details, threads=1):
    #convert format to bam and index
    LOG.write("convertSamToBam(%s, details, %d)\n"%(samFile, threads))
    tempTime = time()
    if 'samtools' not in details['version']:
        command = ["samtools"]
        proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
        proc.wait()
        version_text = proc.stderr.read().decode()
        for line in version_text.splitlines():
            if 'Version' in line:
                details["version"]['samtools'] = line.strip()

    sortThreads = max(int(threads/2), 1)
    samFilePrefix = re.sub(".sam", "", samFile, re.IGNORECASE)
    command = ["samtools", "view", "-bS", "-@", str(sortThreads), "-o", samFilePrefix+"_unsorted.bam", samFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools view return code = %d, time=%d\n"%(return_code, time()-tempTime))
    if return_code != 0:
        comment = "samtools view returned %d"%return_code
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None

    LOG.flush()
    os.remove(samFile) #save a little space

    bamFileSorted = samFilePrefix+".bam" 
    command = ["samtools", "sort", "-@", str(sortThreads), "-o", bamFileSorted, samFilePrefix+"_unsorted.bam"]
    if "Version: 0.1.19" in details["version"]["samtools"]:
        # different invocation for this older version
        command = ["samtools", "sort", "-@", str(sortThreads), samFilePrefix+"_unsorted.bam", samFilePrefix]
    return_code = subprocess.check_call(command, shell=False, stderr=LOG)

    if return_code != 0:
        comment = "samtools sort returned %d, convertSamToBam failed"%return_code
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    LOG.write("bamFileSorted = "+bamFileSorted+"\n")
    if not os.path.exists(bamFileSorted):
        comment = "{0} not found, sorting bamfile failed, convertSamToBam failed\n".format(bamFileSorted)
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    if not os.path.getsize(bamFileSorted):
        comment = "{0} of size zero, sorting bamfile failed, convertSamToBam failed\n".format(bamFileSorted)
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    LOG.write("samtools sort return code=%d, time=%d, size of %s is %d\n"%(return_code, time()-tempTime, bamFileSorted, os.path.getsize(bamFileSorted)))

    command = ["samtools", "index", bamFileSorted]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools index return code = %d\n"%return_code)
    return bamFileSorted

def runRacon(contigFile, longReadsFastq, details, threads=1):
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    Run racon on reads, read-to-contig-sam, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time() - START_TIME))
    LOG.write('runRacon(%s, %s, details, %d)\n'%(contigFile, longReadsFastq, threads))
    if 'racon' not in details['version']:
        command = ["racon", "--version"]
        proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        proc.wait()
        version_text = proc.stderr.read().decode()
        if not version_text:
            version_text = proc.stdout.read().decode()
        details["version"]['racon'] = version_text.strip()
    readsToContigsSam = runMinimap(contigFile, longReadsFastq, details, threads, outformat='sam')
    if not readsToContigsSam:
        comment = "runMinimap failed to generate sam file, exiting runRacon"
        LOG.write(comment + "\n")
        details['problem'].append(comment)
        return None
    raconStartTime = time()
    raconContigs = contigFile.replace(".fasta", ".racon.fasta")
    raconOut = open(raconContigs, 'w')
    command = ["racon", "-t", str(threads), "-u", longReadsFastq, readsToContigsSam, contigFile]
    LOG.write("racon command: \n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL, stdout=raconOut)
    LOG.write("racon return code = %d, time = %d seconds\n"%(return_code, time()-raconStartTime))
    if return_code != 0:
        return None
    raconContigSize = os.path.getsize(raconContigs)
    if raconContigSize < 10:
        return None
    details['polishing'].append({"input_contigs":contigFile, "reads": longReadsFastq, "program": "racon", "output": raconContigs, "seconds": time()-raconStartTime})
    comment = "racon, input %s, output %s"%(contigFile, raconContigs)
    LOG.write(comment+"\n")
    if re.search("racon.racon.fasta", raconContigs):
        shorterFileName = re.sub("racon.racon.fasta", "racon.fasta", raconContigs)
        shutil.move(raconContigs, shorterFileName)
        raconContigs = shorterFileName
        LOG.write("renaming {} to {}\n".format(raconContigs, shorterFileName))
    #details["post-assembly transformation"].append(comment)
    return raconContigs

def runBowtie(contigFile, shortReadFastq, details, threads=1, outformat='bam'):
    """
    index contigsFile, then run bowtie2, then convert sam file to pos-sorted bam and index
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    command = ["bowtie2-build", "--threads", str(threads), contigFile, contigFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout and stderr to dev/null
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("bowtie2-build return code = %d\n"%return_code)
    if return_code != 0:
        return None

    command = ["bowtie2", "-p", str(threads)]
    if details['reads'][shortReadFastq]['platform'] == 'fasta':
        command.append("-f")
    command.extend(["-x", contigFile])
    fastqBase=''
    if ":" in shortReadFastq:
        read1, read2 = shortReadFastq.split(":")
        command.extend(('-1', read1, '-2', read2))
        fastqBase = read1
    else:
        command.extend(('-U', shortReadFastq))
        fastqBase = os.path.basename(shortReadFastq)
    fastqBase = re.sub(r"\..*", "", fastqBase)
    samFile = contigFile+"_"+fastqBase+".sam"
    command.extend(('-S', samFile))
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("bowtie2 return code = %d\n"%return_code)
    if return_code != 0:
        return None
    if outformat == 'sam':
        return samFile

    else:
        contigsBam = convertSamToBam(samFile, details, threads=threads)
        LOG.write('runBowtie returning %s\n'%contigsBam)
        return contigsBam

def runPilon(contigFile, shortReadFastq, details, pilon_jar, threads=1):
    """ 
    polish contigs with short reads (illumina or iontorrent)
    first map reads to contigs with bowtie
    """
    LOG.write("runPilon starting Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    if not (pilon_jar and os.path.exists(pilon_jar)):
        comment = "jarfile %s not found for runPilon, giving up"%(pilon_jar)
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return
    if 'pilon' not in details['version']:
        command = ["java", "-jar", pilon_jar, "--version"]
        proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE)
        proc.wait()
        version_text = proc.stdout.read().decode()
        details["version"]['pilon'] = version_text.strip()

    bamFile = runBowtie(contigFile, shortReadFastq, details, threads=threads, outformat='bam')
    if not bamFile:
        return None
    pilon_start_time = time()
    command = ['java', '-Xmx32G', '-jar', pilon_jar, '--genome', contigFile]
    if ':' in shortReadFastq:
        command.extend(('--frags', bamFile))
    else:
        command.extend(('--unpaired', bamFile))
    pilonPrefix = contigFile.replace(".fasta", "")
    m = re.match(".*pilon_(\d+)", pilonPrefix)
    if m :
        level = int(m.group(1))
        pilonPrefix = re.sub("pilon_"+m.group(1), "pilon_{}".format(level+1), pilonPrefix)
    else:
        pilonPrefix += "_pilon_1"
    pilonContigs = pilonPrefix #+ ".fasta"
    command.extend(('--outdir', '.', '--output', pilonContigs, '--changes'))
    command.extend(('--threads', str(threads)))
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("pilon return code = %d\n"%return_code)
    pilon_time = time() - pilon_start_time
    LOG.write("pilon duration = %d\n"%(pilon_time))
    if return_code != 0:
        return None
    pilon_changes = 0
    with open(pilonContigs+".changes") as CHANGES:
        pilon_changes = len(CHANGES.read().splitlines())
    details['polishing'].append({"input_contigs":contigFile, "reads": shortReadFastq, "program": "pilon", "output": pilonContigs, "num_changes": pilon_changes, "seconds" : pilon_time})
    comment = "pilon, input %s, output %s, num_changes = %d"%(contigFile, pilonContigs, pilon_changes)
    LOG.write(comment+"\n")
    #details["post-assembly transformation"].append(comment)
    return pilonContigs+".fasta" 

def calcReadDepth(bamfiles, details):
    """ Return dict of contig_ids to tuple of (coverage, normalized_coverage) """
    LOG.write("calcReadDepth(%s)\n"%" ".join(bamfiles))
    if 'samtools' not in details['version']:
        command = ["samtools"]
        proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
        proc.wait()
        version_text = proc.stderr.read().decode()
        for line in version_text.splitlines():
            if 'Version' in line:
                details["version"]['samtools'] = line.strip()
    readDepth = {}
    command = ["samtools", "depth"]
    if type(bamfiles) is str:
        command.append(bamfiles)
    else:
        command.extend(bamfiles)
    LOG.write("command = "+" ".join(command)+"\n")
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    depthData = proc.communicate()[0].decode()
    LOG.write("length of depthData string = %d\n"%len(depthData))
    depthSum = 0
    totalDepthSum = 0
    totalLength = 0
    length = 0
    contigLength = {}
    prevContig=None
    for line in iter(depthData.splitlines()):
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            raise Exception("Number of fields is less than 3:\n"+line)
        contig = fields[0]
        depth = 0
        for field in fields[2:]:
            depth += float(field)
        if contig != prevContig:
            if prevContig is not None:
                meanDepth = depthSum/length
                readDepth[prevContig] = [meanDepth, 0]
                contigLength[prevContig] = length
                depthSum = 0
                length = 0
            prevContig = contig
        totalLength += 1
        length += 1
        depthSum += depth
        totalDepthSum += depth
    # process data for last contig
    if prevContig is not None:
        meanDepth = depthSum/length
        readDepth[prevContig] = [meanDepth, 0]
        contigLength[prevContig] = length

    LOG.write("total length for depth data = %d\n"%totalLength)
    LOG.write("total depth = %.1f\n"%totalDepthSum)
    LOG.write("len(readDepth) = %d\n"%len(readDepth))
    totalMeanDepth = 0
    if totalLength > 0:
        totalMeanDepth = totalDepthSum/totalLength
        readDepth['average_coverage'] = totalMeanDepth

    # calculate mean depth of contigs within "normal" boundary around overall mean
    lowerBound = totalMeanDepth * 0.5
    upperBound = totalMeanDepth * 2
    oneXSum = 0.0
    oneXLen = 0
    for c in readDepth:
        meanDepth = readDepth[c][0]
        if meanDepth >= lowerBound and meanDepth <= upperBound:
            oneXSum += meanDepth * contigLength[c]
            oneXLen += contigLength[c]
    oneXDepth = 1
    if oneXLen > 0 and oneXSum > 0:
        oneXDepth = oneXSum/oneXLen # length-weighted average

    for c in readDepth:
        meanDepth = readDepth[c][0]
        normalizedDepth = meanDepth / oneXDepth
        readDepth[c][1] = normalizedDepth
    return readDepth

def runCanu(details, canu_exec="canu", threads=1, genome_size="5m", memory=250, prefix=""):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    canuStartTime = time()
    LOG.write("runCanu: elapsed seconds = %d\n"%(canuStartTime-START_TIME))
    comment = """
usage: canu [-version] [-citation] \
            [-correct | -trim | -assemble | -trim-assemble] \
            [-s <assembly-specifications-file>] \
            -p <assembly-prefix> \
            -d <assembly-directory> \
            genomeSize=<number>[g|m|k] \
            [other-options] \
            [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] file1 file2 ...
    """
    # canu -d /localscratch/allan/canu_assembly -p p6_25X gnuplotTested=true genomeSize=5m useGrid=false -pacbio-raw pacbio_p6_25X.fastq
    # first get canu version
    p = subprocess.Popen([canu_exec, "--version"], shell=False, stdout=subprocess.PIPE)
    canu_version = p.stdout.readline().decode().rstrip()
    details['version']['canu'] = canu_version
    details['assembly']['version'] = canu_version
    details['assembly']['assembler'] = 'canu'
    p.wait()

    command = [canu_exec, "-d", '.', "-p", "canu", "useGrid=false", "genomeSize=%s"%genome_size]
    command.extend(["maxMemory=" + str(memory), "maxThreads=" + str(threads)])
    if "1.7" in canu_version:
        # special handling for this version
        command.append("gnuplotTested=true")
    command.append("stopOnReadQuality=false")
    """
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    pacbio_reads = []
    for item in details['reads']:
        if details['reads'][item]['platform'] in ('pacbio', 'fasta'):
            pacbio_reads.append(item)
            if details['reads'][item]['platform'] == 'fasta':
                comment = 'submitting fasta reads to canu, but calling them "pacbio": '+' '.join(details['platform']['fasta'])
                LOG.write(comment+"\n")
                details['problem'].append(comment)
    if pacbio_reads:
        command.append("-pacbio-raw")
        command.extend(pacbio_reads)
    nanopore_reads = []
    for item in details['reads']:
        if details['reads'][item]['platform'] == 'nanopore':
            nanopore_reads.append(item)
    if nanopore_reads:
        command.append("-nanopore-raw")
        command.extend(nanopore_reads)
    if not pacbio_reads + nanopore_reads:
        LOG.write("no long read files available for canu.\n")
        details["problem"].append("no long read files available for canu")
        return None
    LOG.write("canu command =\n"+" ".join(command)+"\n")
    LOG.flush()

    canuStartTime = time()
    #with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
    with open(os.path.join(DETAILS_DIR, prefix+"canu_stdout.txt"), 'w') as CANU_STDOUT: 
        return_code = subprocess.call(command, shell=False, stdout=CANU_STDOUT, stderr=CANU_STDOUT)
    LOG.write("return code = %d\n"%return_code)
    canuEndTime = time()
    elapsedTime = canuEndTime - canuStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"]['elapsed_time'] = elapsedHumanReadable
    details["assembly"]['command_line'] = " ".join(command)

    LOG.write("Duration of canu run was %s\n"%(elapsedHumanReadable))
    if os.path.exists("canu.report"):
        LOG.write("details_dir = %s\n"%DETAILS_DIR)
        LOG.write("canu_report file name = %s\n"%(prefix+"canu_report.txt"))
        canuReportFile = os.path.join(DETAILS_DIR, (prefix+"canu_report.txt"))
        LOG.write("moving canu.report to %s\n"%canuReportFile)
        shutil.move("canu.report", canuReportFile)
    
    if not os.path.exists("canu.contigs.fasta"):
        comment = "Canu failed to generate contigs file. Check "+prefix+"canu_report.txt"
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        details["problem"].append(comment)
        return None
    # rename to canonical contigs.fasta
    contigsFile = "contigs.fasta"
    shutil.move("canu.contigs.fasta", contigsFile)
    if os.path.exists("canu.contigs.gfa"):
        shutil.move("canu.contigs.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    elif os.path.exists("canu.unitigs.gfa"):
        shutil.move("canu.unitigs.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    details["assembly"]["contigs.fasta size:"] = os.path.getsize(contigsFile)
    return contigsFile

def write_html_report(htmlFile, details):
    LOG.write("writing html report to %s\n"%htmlFile)
    HTML = open(htmlFile, 'w')
    HTML.write("<!DOCTYPE html><html><head>\n")
    HTML.write("<link href=\"https://fonts.googleapis.com/css?family=Work+Sans:300,400,500,600,700,800,900\" rel=\"stylesheet\">\n")
    HTML.write("""
<style>
 body {
    font-family: 'Work Sans', sans-serif;
    color: #444;
    }
 header { padding-bottom: .5in; }
 section { margin-bottom: .5in; }

 a {
    text-decoration: none;
    color: #0d78ef;
    }
 a:hover { text-decoration: underline; }
 h2 { font-size: 1.45em; }
 h2, h3 { font-weight: 500; }
 h2 small { color: #888; font-size: 0.65em; }
 .pull-left { float: left; }
 .pull-right { float: right; }
 sup { display: inline-block; }

 table-num,
 table-ref,
 fig-num,
 fig-ref {
     font-weight: 600;
     }

 /* tables */
 table {
     margin: .1in 0 .5in 0;
     border-collapse: collapse;
     page-break-inside: avoid;
     }
 .table-header {
     color: #fff;
     background: #196E9C;
     padding: 5px;
     text-align: left;
     }
 .table-header th { font-weight: 400; }
 .row-header { border-bottom: 1px solid #196E9C; font-weight: 600; }
 th, td { padding: 5px; }
 th { border-bottom: 2px solid #196E9C; }
 tr:last-child { border-bottom: 1px solid #196E9C; }

 .kv-table,
 .kv-table-2 {
     text-align: left;
     }

 .kv-table-2 td:nth-child(2) {
     border-right: 1px solid rgb(172, 172, 172)
     }

 .kv-table td:first-child,
 .kv-table-2 td:first-child,
 .kv-table-2 td:nth-child(3) {
     font-weight: 700;
     }

 table td.align-right {
     text-align: right;
     }

 .kv-table-2 td:nth-child(2) { padding-right: 10px; }
 .kv-table-2 td:nth-child(3) { padding-left: 10px; }
 .lg-table { width: 80%; }
 .med-table { width: 60%; }
 .sm-table { width: 40%; }
 .xs-table {width: 20%; }

 .center { margin-left: auto; margin-right: auto; }

 .logo {
 width: 2.5in;
 border-right: 3px solid #777;
 padding-right: 10px;
 margin-right: 10px;
 }
 .title {
 padding: 17px 0 0 0px;
 font-size: 1.3em;
 color: #777;
 }
 .report-info {
 text-align: right;
 margin-right: 30px;
 font-size:.8em;
 color: #777;
 }
 .report-info span {
 font-weight: 600;
 }

 .main-img {
 width: 250px;
 margin-right: .2in;
 }

 ol.references li {
 margin-bottom: 1em;
 }

 .clearfix:after {
 visibility: hidden;
 display: block;
 font-size: 0;
 content: " ";
 clear: both;
 height: 0;
 }
 .clearfix { display: inline-block; }

 * html .clearfix { height: 1%; }
 .clearfix { display: block; }
    .a { margin-left: 3em; margin-top: 1em; margin-bottom: 1em}\n")
    HTML.write(".b {margin-left: 3em }\n")
    HTML.write("span.header { font-size: large; font-weight: bold; margin-left: -2em}\n")
 """)
    
    HTML.write("</style></head><body><header>\n")
    HTML.write('<div class="report-info pull-right">\n')
    HTML.write('<span>Report Date:</span> '+strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+'<br>\n')
    HTML.write('</div></header>\n')
    HTML.write("<h1>Genome Assembly Report</h1>\n")

    if "Bandage" in details and "plot" in details["Bandage"]:
        HTML.write('<section>\n<h2>Assembly Plot</h2>\n')
        if os.path.exists(details["Bandage"]["plot"]):

            svg_text = open(details["Bandage"]["plot"]).read()
            svg_text = re.sub(r'<svg width="[\d\.]+mm" height="[\d\.]+mm"', '<svg width="200mm" height="150mm"', svg_text)
            HTML.write("<div class='a'><span class='header'>Bandage Plot:</span><br>\n")
            HTML.write(details["Bandage"]["version"]+"\n")
            HTML.write(svg_text+"\n\n")
            HTML.write("</div>\n")
        HTML.write("</section>\n")

    if 'assembly' in details:
        HTML.write('<section>\n<h2>Assembly</h2>\n')
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Assembly Process </th></tr></thead>
            <tbody>
            """)
        for key in sorted(details['assembly']):
            if key == "problem":
                continue
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['assembly'][key])))
        if "problem" in details['assembly'] and details['assembly']['problem']:
            HTML.write('<tr><td colspan="2">Issues with assembly: </td></tr>')
            for prob in details['assembly']['problem']:
                HTML.write('<tr><td colspan="2">%s:</td></tr>\n'%(prob))
        HTML.write("</tbody></table>\n")
        HTML.write("</section>\n")

    if 'polishing' in details and len(details['polishing']):
        HTML.write('<section>\n<h2>Polishing</h2>\n')
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Polishing Rounds </th></tr></thead>
            <tbody>
            """)
        for iteration, info in enumerate(details['polishing']):
            HTML.write("<tr><td>%s:</td><td>%d</td></tr>\n"%("Round", iteration+1))
            for key in sorted(info):
                HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(info[key])))
            HTML.write("<tr></tr>\n") # blank row
        HTML.write("</tbody></table>\n")
        HTML.write("</section>\n")

    if "contig_filtering" in details and details['contig_filtering']:
        HTML.write('<section>\n<h2>Filtering Contigs on Length and Coverage</h2>\n')
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Contig Filtering </th></tr></thead>
            <tbody>
            """)
        for key in sorted(details['contig_filtering']):
            HTML.write("<tr><td>{}:</td><td>{}</td></tr>\n".format(key, details['contig_filtering'][key])) 
        HTML.write("</tbody></table>\n")
        HTML.write("</section>\n")

    if "quast_txt" in details:
        HTML.write("<section><h2>Quast Report</h2>\n")
        HTML.write("<a href='%s'>%s</a><br>\n"%(details["quast_html"], "Quast html report"))
        HTML.write("</table>\n")
        if os.path.exists(os.path.join(SAVE_DIR, details["quast_txt"])):
            HTML.write("<pre>\n")
            HTML.write(open(os.path.join(SAVE_DIR, details["quast_txt"])).read())
            HTML.write("\n</pre>\n")
        HTML.write("</section>\n")
    
    HTML.write('<section>\n<h2>Input Reads</h2>\n')
    for item in details['reads']:
        if 'supercedes' in details['reads'][item]:
            continue # this is a derived item, not original input
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> """ + item +""" </th></tr></thead>
            <tbody>
            """)
        HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%("read file", item))
        for key in ('platform', 'layout', 'num_reads', 'num_reads_original', 'num_bases', 'num_bases_original',
                'avg_len', 'max_read_len', 'num_read', 'sample_read_id'):  #sorted(details['reads'][item]):
            if key in details['reads'][item]:
                HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
        if "problem" in details['reads'][item] and details['reads'][item]['problem']:
            HTML.write('<tr><td colspan="2">Issues with read set '+item+'</td></tr>')
            for prob in details['reads'][item]['problem']:
                HTML.write('<tr><td colspan="2">%s:</td></tr>\n'%(prob))
    HTML.write("</tbody></table>\n")
    HTML.write("</section>\n")

    if "trim report" in details:
        HTML.write("<section><h2>Trimming Report</h2>\n")
        for reads in details["trim report"]:
            HTML.write("<b>"+reads+"</b><ul>")
            for report in details["trim report"][reads]:
                report_file = os.path.join(DETAILS_DIR, report)
                if os.path.isfile(report_file):
                    HTML.write("<pre>\n")
                    HTML.write(open(report_file).read())
                    HTML.write("\n</pre>\n")
                else:
                    LOG.write("could not open trim report from "+report_file+"\n")
        HTML.write("</section>\n")

    HTML.write("<section><h2>Tools Used:</h2>\n")
    HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Tools Used </th></tr></thead>
            <tbody>
            """)
    for key in sorted(details['version']):
        HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['version'][key])))
    HTML.write("</table>\n")
    HTML.write("</section>\n")
    HTML.write("</html>\n")

    HTML.close()


def main():
    main_return_code = 1 # set to zero when we have an assembly
    global START_TIME
    START_TIME = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='p3_assembly')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    parser.add_argument('--sra', metavar='files', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--max_bases', type=int, default=10000000000, help='process at most this many bases per read file or pair')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('--recipe', choices=['unicycler', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'single-cell', 'rna-spades', 'auto'], help='assembler to use', default='auto')
    parser.add_argument('--contigs', metavar='fasta', help='perform polishing on existing assembly')
    #parser.add_argument('--only-assembler', action='store true', help='omit spades read error correction')
    
    parser.add_argument('--racon_iterations', type=int, default=0, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=0, help='number of times to run pilon per short-read file', required=False)
    parser.add_argument('--pilon_hours', type=float, default=6.0, help='maximum hours to run pilon', required=False)
    #parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--prefix', default='', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', metavar='k, m, or g', default=DEFAULT_GENOME_SIZE, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', type=int, default=300, help='save contigs of this length or longer', required=False)
    parser.add_argument('--min_contig_coverage', type=float, default=5, help='save contigs of this coverage or deeper', required=False)
    #parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=250, help='RAM limit in Gb')
    parser.add_argument('--trim', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--pilon_jar', help='path to pilon executable or jar')
    parser.add_argument('--canu_exec', default="canu", help='path to canu executable (def "canu")')
    parser.add_argument('--spades_for_unicycler', help='path to spades.py suitable for unicycler')
    parser.add_argument('--bandage', action='store_true', help='generate image of assembly path using Bandage')
    parser.add_argument('--params_json', help='JSON file with additional information.')
    parser.add_argument('--path-prefix', help="Add the given directories to the PATH", nargs='*', required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    if args.params_json:
        parseJsonParameters(args)
    baseName = args.outputDirectory #"p3_assembly" 
    if args.prefix:
        args.prefix = args.prefix.replace(" ", "_")
        if not args.prefix.endswith("_"):
            args.prefix += "_"
    global WORK_DIR
    WORK_DIR = baseName+"_work"
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.mkdir(WORK_DIR)
    WORK_DIR = os.path.abspath(WORK_DIR)
    global SAVE_DIR
    SAVE_DIR = os.path.abspath(os.path.join(WORK_DIR, "save"))
    if os.path.exists(SAVE_DIR):
        shutil.rmtree(SAVE_DIR)
    os.mkdir(SAVE_DIR)
    global DETAILS_DIR
    DETAILS_DIR = os.path.abspath(os.path.join(SAVE_DIR, "details"))
    os.mkdir(DETAILS_DIR)
    #logfileName = os.path.join(DETAILS_DIR, args.prefix + "p3_assembly.log")
 
    if args.path_prefix:
        os.environ["PATH"] = ":".join(args.path_prefix) + ":" + os.environ["PATH"]

    global LOG 
    #sys.stderr.write("logging to "+logfileName+"\n")
    #LOG = open(logfileName, 'w') 
    LOG = sys.stderr
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(START_TIME))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Work directory is "+WORK_DIR+"\n\n")
    LOG.write("Final output will be saved to "+SAVE_DIR+"\n\n")
    LOG.write("Detailed output will be saved to "+DETAILS_DIR+"\n\n")
    details = { 'logfile' : logfileName, 'assembly': {} }
    details["pre-assembly transformation"] = []
    details["post-assembly transformation"] = []
    details["original_items"] = []
    details["reads"] = {}
    details["assembly"] = {}
    details["coverage"] = {}
    details["polishing"] = []
    details["version"] = {}
    details["problem"] = []
    details["derived_reads"] = []
    details["platform"] = {'illumina':[], 'iontorrent':[], 'pacbio':[], 'nanopore':[], 'fasta':[], 'anonymous':[]}
    details['max_bases']=args.max_bases

    if args.anonymous_reads:
        categorize_anonymous_read_files(args, details)

    if args.sra:
        fetch_sra_files(args.sra, details)

    if args.illumina:
        platform = 'illumina'
        for item in args.illumina:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved, max_bases=args.max_bases)

    if args.iontorrent:
        platform = 'iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved, max_bases=args.max_bases)

    if args.pacbio:
        for item in args.pacbio:
            registerReads(item, details, platform = 'pacbio', max_bases=args.max_bases)

    if args.nanopore:
        for item in args.nanopore:
            registerReads(item, details, platform = 'nanopore', max_bases=args.max_bases)

    if args.contigs:
        if not os.path.exists(args.contigs):
            comment = "Specified contigs file not found: "+args.contigs
            LOG.write(comment+"\n")
            details['problem'].append(comment)
            sys.stderr.write(comment+"\n")
            exit()
        details['assembly']['input_contigs'] = args.contigs
        comment = "input contigs specified as "+args.contigs
        LOG.write(comment+"\n")
        LOG.write("symlinking %s to %s\n"%(os.path.abspath(args.contigs), os.path.join(WORK_DIR,"input_contigs.fasta")))
        os.symlink(os.path.abspath(args.contigs), os.path.join(WORK_DIR,"input_contigs.fasta"))
    # move into working directory so that all files are local
    os.chdir(WORK_DIR)

    if args.trim and details['platform']['illumina'] + details['platform']['iontorrent']:
        trimGalore(details, threads=args.threads)
    LOG.write("details dir = "+DETAILS_DIR+"\n")
    if args.recipe == "auto":
        #now must decide which assembler to use
        if True:
            # original rule: if any illumina or iontorrent reads present, use Unicycler (long-reads can be present), else use canu for long-reads
            if details['platform']['fasta']:
                args.recipe = "spades"
                LOG.write("auto recipe using spades due to presence of fasta data")
            elif details['platform']['illumina'] + details['platform']['iontorrent']:
                args.recipe = "unicycler"
                LOG.write("auto recipe using unicycler due to presence of short fastq data")
            else:
                args.recipe = "canu"
                LOG.write("auto recipe using canu due to absence of short fastq data")
        else:
            # alternative rule: if any long reads present, use canu
            if details['platform']['pacbio'] + details['platform']['nanopore']:
                args.recipe = "canu"
            else:
                args.recipe = "unicycler"

    contigs = ""
    if args.contigs:
        contigs = "input_contigs.fasta"
            
    elif args.recipe == "unicycler":
        spades_exec = None
        if (args.spades_for_unicycler):
            spades_exec = args.spades_for_unicycler
        contigs = runUnicycler(details, threads=args.threads, min_contig_length=args.min_contig_length, prefix=args.prefix, spades_exec=spades_exec )
        if not contigs:
            comment = "unicycler failed to generate contigs, trying spades"
            LOG.write(comment+"\n")
            details['problem'].append(comment)
            #details["pre-assembly transformation"].append(comment)
            contigs = runSpades(details, args)
    elif args.recipe == "canu":
        contigs = runCanu(details, canu_exec=args.canu_exec, threads=args.threads, genome_size=args.genome_size, memory=args.memory, prefix=args.prefix)
    elif "spades" in args.recipe or args.recipe == "single-cell":
        if args.recipe == "meta-spades" and (args.pilon_iterations or args.racon_iterations):
            args.pilon_iterations = 0
            args.racon_iterations = 0
            comment = "Because recipe is meta-spaces, turning pilon and racon iterations off."
            LOG.write(comment+"\n")
            details['problem'].append(comment)
        contigs = runSpades(details, args)
    else:
        LOG.write("cannot interpret args.recipe: "+args.recipe)

    if contigs and os.path.getsize(contigs):
        if not args.contigs:
            main_return_code = 0
        LOG.write("size of contigs file is %d\n"%os.path.getsize(contigs))
        if args.racon_iterations:
            # now run racon with each long-read file
            for longReadFile in details['reads']:
                if details['reads'][longReadFile]['platform'] == 'fasta':
                    continue # do not run racon on fasta reads
                if details['reads'][longReadFile]['length_class'] == 'long':
                    try:
                        for i in range(0, args.racon_iterations):
                            LOG.write("runRacon(%s, %s, round=%d, details, threads=%d)\n"%(contigs, longReadFile, i, args.threads))
                            raconContigFile = runRacon(contigs, longReadFile, details, threads=args.threads)
                            if raconContigFile is not None:
                                contigs = raconContigFile
                            else:
                                break # break out of iterating racon_iterations, go to next long-read file if any
                    except Exception as e:
                        comment = "runRacon failed with exception {}".format(e)
                        LOG.write(comment)
                        sys.stderr.write(comment)
            
        if args.pilon_iterations and args.pilon_jar:

            pilon_end_time = time() + args.pilon_hours * 60 * 60
            # now run pilon with each short-read file
            for readFastq in details['reads']:
                if time() > pilon_end_time:
                    break
                if 'superceded_by' in details['reads'][readFastq]:
                    continue # may have been superceded by trimmed version of those reads
                if details['reads'][readFastq]['platform'] == 'fasta':
                    continue # do not run pilon on fasta reads
                if details['reads'][readFastq]['length_class'] == 'short':
                    try:
                        for iteration in range(0, args.pilon_iterations):
                            if time() > pilon_end_time:
                                LOG.write("Time expended on pilon exceeds allocation of {} hours. Omitting further pilon runs.".format(args.pilon_hours))
                                break

                            LOG.write("runPilon(%s, %s, round=%d, details, threads=%d)\n"%(contigs, readFastq, iteration, args.threads))
                            pilonContigFile = runPilon(contigs, readFastq, details, args.pilon_jar, args.threads)
                            if pilonContigFile is not None:
                                contigs = pilonContigFile
                            else:
                                break
                        # check number of changes in most recent round, quit if zero
                        if 'num_changes' in details['polishing'][-1] and details['polishing'][-1]['num_changes'] == 0:
                            break
                    except Exception as e:
                        comment = "runPilon failed with exception {}".format(e)
                        LOG.write(comment)
                        sys.stderr.write(comment)
            
    if contigs and os.path.getsize(contigs):
        filteredContigs = filterContigsByMinLength(contigs, details, args.min_contig_length, args.min_contig_coverage, args.threads, args.prefix)
        if filteredContigs:
            contigs = filteredContigs
    if contigs and os.path.getsize(contigs):
        runQuast(contigs, args, details)
        shutil.move(contigs, os.path.join(SAVE_DIR, args.prefix+"contigs.fasta"))

    gfaFile = os.path.join(DETAILS_DIR, args.prefix+"assembly_graph.gfa")
    if os.path.exists(gfaFile) and os.path.getsize(gfaFile):
        runBandage(gfaFile, details)

    with open(os.path.join(DETAILS_DIR, args.prefix+"run_details.json"), "w") as fp:
        try:
            json.dump(details, fp, indent=2, sort_keys=True)
        except UnicodeDecodeError as ude:
            LOG.write("Problem writing details to json: "+str(ude)+"\n")


    htmlFile = os.path.join(SAVE_DIR, args.prefix+"assembly_report.html")
    write_html_report(htmlFile, details)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - START_TIME)/3600))
    LOG.write("Total time in seconds = %d\t"%((time() - START_TIME)))
    LOG.close()
    return main_return_code

if __name__ == "__main__":
    main()
