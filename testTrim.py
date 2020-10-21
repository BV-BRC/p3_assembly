#!/usr/bin/env python
import sys
import subprocess
import argparse
from tempfile import mkdtemp
import gzip
import os
import os.path
import re
import shutil
import urllib2
from time import time, localtime, strftime
import json
#import sra_tools
import glob
import StringIO


def registerReads(reads, details, platform=None, interleaved=False, supercedes=None):
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
    
    readStruct = {}
    readStruct["file"] = []
    readStruct["path"] = []
    readStruct["problems"] = []
    readStruct["layout"] = 'na'
    readStruct["platform"] = 'na'
    if ":" in reads or "%" in reads:
        if ":" in reads:
            delim = ["%", ":"][":" in reads] # ":" if it is, else "%"
        read1, read2 = reads.split(delim)
        readStruct["delim"] = delim
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
        if os.path.abspath(dir1) != WorkDir:
            LOG.write("symlinking %s to %s\n"%(read1, os.path.join(WorkDir,file1)))
            os.symlink(os.path.abspath(read1), os.path.join(WorkDir,file1))
        if os.path.abspath(dir2) != WorkDir:
            LOG.write("symlinking %s to %s\n"%(read2, os.path.join(WorkDir,file2)))
            os.symlink(os.path.abspath(read2), os.path.join(WorkDir,file2))
        readStruct["file"].append(file1)
        readStruct["file"].append(file2)
        readStruct["path"].append(read1)
        readStruct["path"].append(read2)
        # the "files" entry is an array1 of tuples, each tuple is a path, basename yielded by os.path.split
        registeredName = delim.join(sorted((file1, file2)))
    else:
        # no ':' or '%' delimiter, so a single file
        if not os.path.exists(reads):
            comment = "file does not exist: %s"%reads
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        if interleaved:
            readStruct["interleaved"] = True
        dir1, file1 = os.path.split(reads)
        if os.path.abspath(dir1) != WorkDir:
            LOG.write("symlinking %s to %s\n"%(reads, os.path.join(WorkDir,file1)))
            os.symlink(os.path.abspath(reads), os.path.join(WorkDir,file1))
        readStruct["file"].append(file1)
        readStruct["path"].append(reads)
        registeredName = file1

    if supercedes:
        details['derived_reads'].append(registeredName)
        readStruct['supercedes'] = supercedes
        if supercedes in details['reads']:
            details['reads'][supercedes]['superceded_by'] = registeredName
            platform = details['reads'][supercedes]['platform']
            readStruct['platform'] = platform
            try: # swap item in platform[] with this new version
                index = details['platform'][platform].index(supercedes)
                details['platform'][platform][index] = registeredName
            except ValueError as ve:
                comment = "Problem: superceded name %s not found in details_%s"%(registeredName, platform)
    else:
        if platform:
            readStruct['platform'] = platform
            if platform not in details['platform']:
                details["platform"][platform] = []
            details["platform"][platform].append(registeredName)
    if registeredName in details['reads']:
        comment = "registered name %s already in details[reads]"%registeredName
        LOG.write(comment+"\n")
        details["problem"].append(comment)

    details['reads'][registeredName]=readStruct
    if len(readStruct['file']) == 2:
        studyPairedReads(registeredName, details)
    else:
        studySingleReads(registeredName, details)
    return registeredName

def parseJasonParameters(args):
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
        if "srr_ids" in data:
            if not args.sra:
                pass
    return

def studyPairedReads(item, details):
    """
    Read both files. Verify read ID are paired. Determine avg read length. Update details['reads'].
    """
    func_start = time()
    LOG.write("studyPairedReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start - Start_time))
    if "reads" not in details:
        details["reads"] = {}
    if item not in details["reads"]:
        details["reads"][item] = {}
    details['reads'][item]['layout'] = 'paired-end'
    details['reads'][item]['avg_len'] = 0
    details['reads'][item]['length_class'] = 'na'
    details['reads'][item]['num_reads'] = 0
    file1, file2 = item.split(":")
    if file1.endswith("gz"):
        F1 = gzip.open(os.path.join(WorkDir, file1))
        F2 = gzip.open(os.path.join(WorkDir, file2))
    else:
        F1 = open(os.path.join(WorkDir, file1))
        F2 = open(os.path.join(WorkDir, file2))

    line = F1.readline()
    sample_read_id = line.split(' ')[0]
    F1.seek(0)
    if sample_read_id.startswith('>'):
        details['reads'][item]['platform']='fasta'
        if 'fasta' not in details['platform']:
            details['platform']['fasta'] = []
        details['platform']['fasta'].append(file1)
        details['platform']['fasta'].append(file2)
        studyFastaReads(F1, file1, details)
        studyFastaReads(F2, file2, details)
        F1.close()
        F2.close()
        return
    
    read_ids_paired = True
    seqLen1 = 0
    seqLen2 = 0
    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    i = 0
    for line1 in F1:
        line2 = F2.readline()
        if not line2:
            line2 = ""
        if i % 4 == 0 and read_ids_paired:
            read_id_1 = line1.split(' ')[0] # get part up to first space, if any 
            read_id_2 = line2.split(' ')[0] # get part up to first space, if any 
            if not readNumber:
                sample_read_id = read_id_1
            if not read_id_1 == read_id_2:
                diff = findSingleDifference(read_id_1, read_id_2)
                if diff == None or sorted(diff) != ('1', '2'):
                    read_ids_paired = False
                    details['reads'][item]["problems"].append("id_mismatch at read %d: %s vs %s"%(readNumber+1, read_id_1, read_id_2))
        elif i % 4 == 1:
            seqLen1 = len(line1)-1
            seqLen2 = len(line2)-1
        elif i % 4 == 3:
            if seqQualLenMatch:
                if not (seqLen1 == len(line1)-1 and seqLen2 == len(line2)-1):
                    readId = [read_id_1, read_id_2][seqLen1 != len(line1)-1]
                    seqQualLenMatch = False
                    comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, readId)
                    details['reads'][item]["problems"].append(comment)
                    LOG.write(comment+"\n")
            totalReadLength += seqLen1 + seqLen2
            maxReadLength = max(maxReadLength, seqLen1, seqLen2) 
            minReadLength = min(minReadLength, seqLen1, seqLen2)
            minQualScore = min(minQualScore + line1.rstrip() + line2.rstrip())
            maxQualScore = max(maxQualScore + line1.rstrip() + line2.rstrip())
            readNumber += 1
        i += 1

    F1.close()
    F2.close()

    avgReadLength = totalReadLength/(readNumber*2)
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 

    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, details['reads'][item]['avg_len'])
    details['reads'][item]['length_class'] = ["short", "long"][avgReadLength >= Max_short_read_length]
    if avgReadLength >= Max_short_read_length:
        comment = "paired reads appear to be long, expected short: %s"%item
        LOG.write(comment+"\n")
        details['reads'][item]['problem'].append(comment)

    LOG.write("duration of studyPairedReads was %d seconds\n"%(time() - func_start))
    return

def studySingleReads(item, details):
    func_start = time()
    LOG.write("studySingleReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start-Start_time))
    if "reads" not in details:
        details["reads"] = {}
    if item not in details['reads']:
        details["reads"][item] = {}
    details['reads'][item]['layout'] = 'single-end'
    details['reads'][item]['num_reads'] = 0
    details['reads'][item]['problems'] = []
    if item.endswith("gz"):
        F = gzip.open(item)
    else:
        F = open(item)

    line = F.readline()
    sample_read_id = line.split(' ')[0] # get part up to first space, if any 
    F.seek(0)
    if sample_read_id.startswith(">"):
        studyFastaReads(F)
        if 'fasta' not in details['platform']:
            details['platform']['fasta'] = []
        details['platform']['fasta'].append(item)
        details['reads'][item]['length_class'] = 'long'
        details['reads'][item]['platform'] = 'fasta'
        studyFastaReads(F1, item, details)
        F.close()
        return

    seqLen1 = 0
    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    i = 0
    interleaved = True
    prev_read_id = None
    for line in F:
        if i % 4 == 0:
            read_id = line.split(' ')[0] # get part up to first space, if any 
            if not sample_read_id:
                sample_read_id = read_id
            if interleaved and i % 8 == 0 and prev_read_id: # at every other sample ID check for matching prev, indicates interleaved
                if prev_read_id != sample_read_id:
                    diff = findSingleDifference(prev_read_id, sample_read_id)
                    if diff == None or sorted(diff) != ('1', '2'):
                        interleaved=False

        elif i % 4 == 1:
            seqLen = len(line)-1
        elif i % 4 == 3:
            qualLen = len(line)-1
            if seqQualLenMatch and (seqLen != qualLen):
                seqQualLenMatch = False
                comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, read_id)
                details['reads'][item]["problems"].append(comment)
                LOG.write(comment+"\n")
            totalReadLength += seqLen
            maxReadLength = max(maxReadLength, seqLen) 
            minReadLength = min(minReadLength, seqLen)
            minQualScore = min(minQualScore + line.rstrip())
            maxQualScore = max(maxQualScore + line.rstrip())
            readNumber += 1
        i += 1
                
    avgReadLength = totalReadLength/readNumber
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 
    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, avgReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][avgReadLength >= Max_short_read_length]
    if interleaved:
        details['reads'][item]['interleaved'] = True

    LOG.write("duration of studySingleReads was %d seconds\n"%(time() - func_start))
    return

def studyFastaReads(F, item, details):
    """
    assume format is fasta
    F is opened file
    count reads, calc total length, mean, max, min
    """
    seq = ""
    seqLen = 0
    totalReadLength = 0
    maxReadLength = 0
    minReadLength = 1e6
    readNumber = 0
    i = 0
    for line in F:
        if line.startswith(">"):
            if seq:
                seqLen = len(seq)
                totalReadLength += seqLen
                maxReadLength = max(maxReadLength, seqLen) 
                minReadLength = min(minReadLength, seqLen)
                readNumber += 1
                seq = ""
        else:
            seq += line.rstrip()
    if seq:
        seqLen = len(seq)
        totalReadLength += seqLen
        maxReadLength = max(maxReadLength, seqLen) 
        minReadLength = min(minReadLength, seqLen)
        readNumber += 1

    avgReadLength = totalReadLength/readNumber
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 
    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, avgReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][avgReadLength >= Max_short_read_length]

    LOG.write("duration of studyFastaReads was %d seconds\n"%(time() - func_start))
    return


def inferPlatform(read_id, avgReadLength):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    if read_id.startswith(">"):
        return "fasta"
    if avgReadLength < Max_short_read_length:
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
# NOTE: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
    if re.match(r"@\S+/\S+/\S+_\S+$", read_id): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
        return "pacbio" # 
    if re.match(r"@\S+/\S+$", read_id): #@<MovieName>/<ZMW_number> 
        return "pacbio" # 
#@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
    if re.match(r"@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", read_id): #based on one example, need to test more 
        return "nanopore" # 
    return "pacbio" # default long fastq type

def trimGalore(details, threads=1):
    startTrimTime = time()
    LOG.write("\ntrimGalore() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    if "trim report" not in details:
        details["trim report"] = {}
    toRegister = {} # save trimmed reads to register after iterating dictionary to avoid error "dictionary changed size during iteration"
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'short' and not details['reads'][reads]['platform'] == 'fasta':
            command = ['trim_galore', '-j', str(threads), '-o', '.']
            if ':' in reads:
                splitReads = reads.split(":")
                #-j 4 -o testTrim --length 30 --paired SRR1395326_1_10pct.fastq SRR1395326_2_10pct.fastq
                command.extend(["--paired", splitReads[0], splitReads[1]])
                LOG.write("command: "+" ".join(command)+"\n")
                proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
                trimGaloreStderr = proc.stderr.read()
                return_code = proc.wait()
                LOG.write("return code = %d\n"%return_code)
                trimReads = re.findall(r"Writing validated paired-end read \d reads to (\S+)", trimGaloreStderr)
                LOG.write("re.findall for trimmed files returned %s\n"%str(m))
                if not trimReads or len(trimReads) < 2:
                    comment = "trim_galore did not name trimmed reads output files in stderr"
                    LOG.write(comment+"\n")
                    details['reads'][reads]['problem'].append(comment)
                    continue
                comment = "trim_galore, input %s, output %s"%(reads, ":".join(trimReads))
                LOG.write(comment+"\n")
                details["pre-assembly transformation"].append(comment)
                toRegister[trimmedReadPair] = trimReads[0]+":"+trimReads[1]

                trimReports = re.findall(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
                LOG.write("re.findall for trim reports returned %s\n"%str(m))
                details["trim report"][reads]=[]
                for f in trimReports:
                    shutil.move(f, os.path.join(SaveDir, os.path.basename(f)))
                    details["trim report"][reads].append(f)
            else:
                command.append(reads)
                LOG.write("command: "+" ".join(command))
                return_code = subprocess.call(command, shell=False, stderr=LOG)
                LOG.write("return code = %d\n"%return_code)
                trimReads = reads
                trimReport = trimReads+"_trimming_report.txt"
                trimReads = re.sub(r"(.*)\..*", r"\1_trimmed.fq", trimReads)
                comment = "trim_galore, input %s, output %s"%(reads, trimReads)
                if os.path.exists(trimReads):
                    toRegister[trimReads] = reads
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                    if os.path.exists(trimReport):
                        shutil.move(trimReport, os.path.join(SaveDir, trimReport))
                        details["trim report"][reads]=[trimReport]
                else:
                    comment = "Problem during trim_galore: expected files not found: "+trimReads
                    LOG.write(comment)
                    details["problem"].append(comment)
    for trimReads in toRegister:
        registerReads(trimReads, details, supercedes=toRegister[trimReads])

if __name__ == "__main__":
    details = {}
    details["pre-assembly transformation"] = []
    details["post-assembly transformation"] = []
    details["original_items"] = []
    details["reads"] = {}
    details["problem"] = []
    details["derived_reads"] = []
    details["platform"] = {'illumina':[], 'iontorrent':[], 'pacbio':[], 'nanopore':[], 'fasta':[], 'anonymous':[]}
    global LOG
    global WorkDir
    WorkDir = os.path.abspath(".")
    global Start_time
    Start_time = time()
    global Max_short_read_length
    Max_short_read_length = 400
    LOG = open(WorkDir +"/testTrim.log", "w")

    reads = sys.argv[1]
    reads = registerReads(reads, details, platform='illumina')
    LOG.write("before trimGalore, details = \n%s\n\n"%str(details))
    trimGalore(details)
