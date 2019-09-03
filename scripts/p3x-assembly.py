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

Default_genome_size = "5m"
Default_bytes_to_sample = 20000
Default_window_width = 4
Default_window_quality = 15
Default_end_quality = 10
Max_short_read_length = 600
Read_id_sample = {}
Read_file_type = {}
Avg_read_length = {}
LOG = None # create a log file at start of main()
Start_time = None
WorkDir = None
SaveDir = None

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
                
    if not readNumber:
        comment = "no reads found in %s\n"%item
        LOG.write(comment+"\n")
        return
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
                LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
                if not trimReads or len(trimReads) < 2:
                    comment = "trim_galore did not name trimmed reads output files in stderr"
                    LOG.write(comment+"\n")
                    details['reads'][reads]['problem'].append(comment)
                    continue
                comment = "trim_galore, input %s, output %s"%(reads, ":".join(trimReads))
                LOG.write(comment+"\n")
                details["pre-assembly transformation"].append(comment)
                toRegister[":".join(trimReads)] = reads

                trimReports = re.findall(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
                LOG.write("re.findall for trim reports returned %s\n"%str(trimReports))
                details["trim report"][reads]=[]
                for f in trimReports:
                    shutil.move(f, os.path.join(SaveDir, os.path.basename(f)))
                    details["trim report"][reads].append(f)
            else:
                command.append(reads)
                LOG.write("command: "+" ".join(command))
                proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
                trimGaloreStderr = proc.stderr.read()
                return_code = proc.wait()
                LOG.write("return code = %d\n"%return_code)
                trimReads = re.search(r"Writing final adapter and quality trimmed output to (\S+)", trimGaloreStderr)
                LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
                if not trimReads:
                    comment = "trim_galore did not name trimmed reads output files in stderr"
                    LOG.write(comment+"\n")
                    details['reads'][reads]['problem'].append(comment)
                    continue
                trimReads = trimReads.group(1)
                comment = "trim_galore, input %s, output %s"%(reads, trimReads)
                LOG.write(comment+"\n")
                details["pre-assembly transformation"].append(comment)
                toRegister[trimReads] = reads

                trimReport = re.search(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
                LOG.write("regex for trim report returned %s\n"%str(trimReport))
                if trimReport:
                    trimReport = trimReport.group(1)
                    details["trim report"][reads]=trimReport
                    shutil.move(trimReport, os.path.join(SaveDir, os.path.basename(trimReport)))
    for trimReads in toRegister:
        registerReads(trimReads, details, supercedes=toRegister[trimReads])

    LOG.write("trim_galore trimming completed, duration = %d seconds\n\n\n"%(time()-startTrimTime))

def sampleReads(filename, details=None):
    srf_time = time()
    LOG.write("sampleReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(srf_time)), srf_time-Start_time))
    # figures out Read_file_type
    #return read_format and sample of read ids
    read_format = 'na'
    read_id_sample = []

    if filename.endswith("gz"):
        F = gzip.open(filename)
    else:
        F = open(filename)
    text = F.read(Default_bytes_to_sample) #read X number of bytes for text sample
    F.close()

    LOG.write("  file text sample %s:\n%s\n\n"%(filename, text[0:50]))
    lines = text.split("\n")
    if len(lines) < 2:
        comment = "in sampleReads for %s: text sample (length %d) lacks at least 2 lines"%(filename, len(text))
        LOG.write(comment+"\n")
        details.problem.append(comment)
    if lines[0].startswith("@"):
        read_format = 'fastq'
    elif lines[0].startswith(">"):
        read_format = 'fasta'
        read_id_sample.append(lines[0].split()[0])
    if read_format == 'fastq':
        avg_read_length = 0
        readLengths = []
        for i, line in enumerate(lines):
            if i % 4 == 0:
                read_id_sample.append(line.split(' ')[0]) # get part up to first space, if any 
            elif i % 4 == 1:
                readLengths.append(len(line)-1)
        avg_read_length = sum(readLengths)/float(len(readLengths))
        LOG.write("fastq read type %s average read length %.1f\n"%(read_format, avg_read_length))
    return read_id_sample, avg_read_length

def findSingleDifference(s1, s2):
# if two strings differ in only a single position, return the chars at that pos, else return None
    if len(s1) != len(s2):
        return None
    retval = None
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            if retval:
                return None
            retval = (c1, c2)
    return retval

def categorize_anonymous_read_files(args, details):
    LOG.write("categorize_anonymous_read_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("  files=%s\n"%("\t".join(args.anonymous_reads)))
    read_file_type = {}
    read_id_sample = {}
    singleFiles = []
    pairedFiles = []
    for item in args.anonymous_reads:
        if ":" in item or "%" in item:
            pairedFiles.append(item)
        else:
            read_id_sample[item], avg_read_length = sampleReads(item, details)
            read_file_type[item] = inferPlatform(read_id_sample[item][0], avg_read_length)
            comment = "interpreting %s type as %s"%(item, read_file_type[item])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)
            if read_file_type[item] is not None:
                singleFiles.append(item)

    # try to find paired files
    membersOfPairs = set()
    for i, filename1 in enumerate(singleFiles[:-1]):
        for filename2 in singleFiles[i+1:]:
            singleDiff = findSingleDifference(filename1, filename2)
# singleDiff will be not None if the strings match at all but one character (presumably '1' vs '2')
            if singleDiff:
                intDiffs = (int(singleDiff[0]), int(singleDiff[1]))
                pair = None
                if intDiffs[0] == 1 and intDiffs[1] == 2:
                    pair = (filename1, filename2)
                elif intDiffs[1] == 1 and intDiffs[0] == 2:
                    pair = (filename2, filename1)
                if pair:
                    comment = "candidate paired files: %s  %s"%pair
                    LOG.write(comment+"\n")
                    details.problem.append(comment)
                    if read_file_type[filename1] != read_file_type[filename2]:
                        comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
                        LOG.write(comment+"\n")
                        details.problem.append(comment)
                        continue
                    pairedFiles.append(pair)

    # now go over all pairs to test for matching types and matching read IDs
    valid_pairs = set()
    valid_singles = set()
    for item in pairedFiles:
        if ":" in item:
            filename1, filename2 = item.split(":") 
        elif "%" in item:
            filename1, filename2 = item.split("%") 
        else:
            comment = "failed to find ':' or '%' in file pair: "+item
            LOG.write(comment+"\n")
            details.problem.append(comment)
            continue
        if filename1 not in read_file_type:
            read_id_sample[filename1], avg_read_length = sampleReads(filename1, details)
            read_file_type[filename1] = inferPlatform(read_id_sample[filename1][0], avg_read_length)
            comment = "interpreting %s type as %s"%(filename1, read_file_type[filename1])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)
        if filename2 not in read_file_type:
            read_id_sample[filename2], avg_read_length = sampleReads(filename2, details)
            read_file_type[filename2] = inferPlatform(read_id_sample[filename2][0], avg_read_length)
            comment = "interpreting %s type as %s"%(filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)

        read_types_match = True
        # test if read types are the same
        if read_file_type[filename1] != read_file_type[filename2]:
            comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            details.problem.append(comment)
            read_types_match = False
        read_file_type[item] = read_file_type[filename1] #easier to retrieve later

        # test if read IDs match between files
        ids_paired = True
        for idpair in zip(read_id_sample[filename1], read_id_sample[filename2]):
            if idpair[0] != idpair[1] and not findSingleCharDiff(idpair[0], idpair[1]):
                ids_paired = False
                comment = "Read IDs do not match for %s(%s) vs %s(%s)"%(filename1, idpair[0], filename2, idpair[1])
                LOG.write(comment+"\n")
                details.problem.append(comment)
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
        registerReads(item, details, platform=read_file_type[item], interleaved = args.interleaved and item in args.interleaved)

    return

def get_sra_runinfo(run_accession, log=None):
    """ take sra run accession (like SRR123456)
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    Altered from versionin sra_tools to handle case of multiple sra runs returned by query.
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

def fetch_one_sra(sra, run_info=None, log=sys.stderr):
    """ requires run_info to know which program to use
    """
    if not run_info:
        run_info = get_sra_runinfo(sra, log)
    programToUse = "fasterq-dump" # but not appropriate for pacbio or nanopore
    if run_info['Platform'].startswith("PACBIO") or run_info['Platform'].startswith("OXFORD_NANOPORE"):
        programToUse = 'fastq-dump'
    stime = time()
    command = [programToUse, "--split-files", sra]
    log.write("command = "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=log)
    log.write("return_code = %d\n"%(return_code))
    if return_code != 0:
        log.write("Problem, %s return code was %d\n"%(programToUse, return_code))

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
    LOG.write("fetch_sra_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
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
        item = None
        if runinfo['LibraryLayout'].startswith("PAIRED"):
            if len(fastqFiles) == 2:
                item = ":".join(sorted(fastqFiles)[:2])
            else:
                comment = "for PAIRED library %s, number of files was %s, expected 2"%(sra, len(fastqFiles))
                details['problem'].append(comment)
                LOG.write(comment+"\n")
        if not item: # library layout single or failed above
            if len(fastqFiles) == 1:
                item = fastqFiles[0]
            if len(fastqFiles) > 1:
                subprocess.call("cat %s*fastq > %s.fastq"%(sra, sra), shell=True)
                item = sra+".fastq"
                comment = "for library %s, list of files was %s, concatenated to %s"%(sra, str(fastqFiles), item)
                details['problem'].append(comment)
                LOG.write(comment+"\n")
            if not item:
                comment = "for %s no fastq file found"%sra
                details['problem'].append(comment)
                LOG.write(comment+"\n")
                continue # failed on that sra
   
        platform = None
        if runinfo["Platform"] == "ILLUMINA":
            platform = "illumina"
        elif runinfo["Platform"] == "ION_TORRENT":
            platform = "iontorrent"
        elif runinfo["Platform"] == "PACBIO_SMRT":
            platform = "pacbio"
        elif runinfo["Platform"] == "OXFORD_NANOPORE":
            platform = "nanopore"
        if not platform:
            ids, avg_length = sampleReads(fastqFiles[0], details)
            platform = ["illumina", "pacbio"][avg_length > Max_short_read_length]
        registerReads(item, details, platform=platform)

    return

def writeSpadesYamlFile(details):
    LOG.write("writeSpadesYamlFile: elapsed seconds = %f\n"%(time()-Start_time))
    outfileName = "spades_yaml_file.txt"
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    for platform in details['platform']:
        LOG.write(platform+": "+", ".join(details['platform'][platform])+"\n")
    
    single_end_reads = []
    paired_end_reads = [[], []]
    mate_pair_reads = [[], []]
    interleaved_reads = []
    all_read_files = []

    shortReadItems = []
    if 'illumina' in details['platform']:
        shortReadItems.extend(details['platform']['illumina'])
    if 'iontorrent' in details['platform']:
        shortReadItems.extend(details['platform']['iontorrent'])
    for item in shortReadItems:
        if ":" in item:
            f = details['reads'][item]['file'][0]
            paired_end_reads[0].append(f)
            f = details['reads'][item]['file'][1]
            paired_end_reads[1].append(f)
        elif "%" in item:
            f = details['reads'][item]['file'][0]
            mate_pair_reads[0].append(f)
            f = details['reads'][item]['file'][1]
            mate_pair_reads[1].append(f)
        else:
            f = details['reads'][item]['file'][0]
            if 'interleaved' in details['reads'][item]:
                interleaved_reads.append(f)
            else:
                single_end_reads.append(f)

    precedingElement=False
    if len(single_end_reads):
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(single_end_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if len(interleaved_reads):
        OUT.write("  {\n    type: \"paired-end\",\n    interlaced reads: [\n        \"")
        OUT.write("\",\n        \"".join(interleaved_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if len(paired_end_reads[0]):
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"fr\",\n")
        OUT.write("    type: \"paired-end\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(paired_end_reads[0]))
        OUT.write("\"\n    ],\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(paired_end_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement=True
    if len(mate_pair_reads[0]):
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"rf\",\n")
        OUT.write("    type: \"mate-pairs\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[0]))
        OUT.write("\"\n    ]\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement=True
    if len(details['platform']['pacbio']):
        pacbio_reads = []
        for item in details['platform']['pacbio']:
            f = details['reads'][item]['file'][0]
            pacbio_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if len(details['platform']['nanopore']):
        nanopore_reads = []
        for item in details['platform']['nanopore']:
            f = details['reads'][item]['path'][0]
            nanopore_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if len(details['platform']['fasta']):
        fasta_reads = []
        for item in details['platform']['fasta']:
            f = details['reads'][item]['file'][0]
            fasta_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"untrusted-contigs\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(fasta_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True

    OUT.write("]\n")
    OUT.close()
    return(outfileName)    

def runQuast(contigsFile, args, details):
    LOG.write("runQuast() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    quastDir = "quast_out"
    quastCommand = ["quast.py",
                    "-o", quastDir,
                    "-t", str(args.threads),
                    "--min-contig", str(args.min_contig_length),
                    contigsFile]
    LOG.write("running quast: "+" ".join(quastCommand)+"\n")
    return_code = subprocess.call(quastCommand, shell=False, stderr=LOG)
    LOG.write("return code = %d\n"%return_code)
    shutil.move(os.path.join(quastDir, "report.html"), os.path.join(SaveDir, args.prefix+"quast_report.html"))
    shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(SaveDir, args.prefix+"quast_report.txt"))
    #shutil.move(os.path.join(quastDir, "icarus_viewers"), os.path.join(SaveDir, "icarus_viewers"))
    details["quast_txt"] = args.prefix+"quast_report.txt"
    details["quast_html"] = args.prefix+"quast_report.html"

def filterContigsByMinLength(inputFile, args, details, shortReadDepth=None, longReadDepth=None):
    """ 
    Write only sequences at or above min_length to output file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("filterContigsByMinLength(%s) \n"%(inputFile))
    if not shortReadDepth and not longReadDepth:
        LOG.write("Neither shortReadDepth nor longReadDepth provided. Quiting.\n")
        return None
    suboptimalReads = ""
    outputFile = re.sub(r"\..*", "_depth_cov_filtered.fasta", inputFile)
    LOG.write("writing filtered contigs to %s\n"%outputFile)
    with open(inputFile) as IN:
        with open(outputFile, 'w') as OUT:
            shortCovTags = ['coverage', 'normalized_cov']
            longCovTags = ['coverage_longread', 'normalized_cov_longread']
            seqId=None
            seq = ""
            contigIndex = 1
            contigCoverage = 0
            contigInfo = ""
            for line in IN:
                m = re.match(r">(\S+)", line)
                if m:
                    if seq:
                        contigId = ">"+args.prefix+"contig_%d length %5d"%(contigIndex, len(seq))
                        contigIndex += 1
                        contigCoverage = 0
                        if shortReadDepth and seqId in shortReadDepth:
                            meanDepth, normalizedDepth = shortReadDepth[seqId]
                            contigId += " %s %.01f %s %.2f"%(shortCovTags[0], meanDepth, shortCovTags[1], normalizedDepth)
                            contigCoverage = meanDepth
                        if longReadDepth and seqId in longReadDepth:
                            meanDepth, normalizedDepth = longReadDepth[seqId]
                            contigId += " %s %.01f %s %.2f"%(longCovTags[0], meanDepth, longCovTags[1], normalizedDepth)
                            contigCoverage = max(meanDepth, contigCoverage)
                        if "contigCircular" in details and contigIndex in details["contigCircular"]:
                            contigId += " circular=true"
                        if len(seq) >= args.min_contig_length and contigCoverage >= args.min_contig_coverage:
                            OUT.write(contigId+"\n")
                            for i in range(0, len(seq), 60):
                                OUT.write(seq[i:i+60]+"\n")
                        else:
                            suboptimalReads += contigId+"\n"+seq+"\n"
                    seq = ""
                    seqId = m.group(1)
                else:
                    seq += line.rstrip()
            if seq:
                contigId = ">"+args.prefix+"contig_%d length %5d"%(contigIndex, len(seq))
                contigCoverage = 0
                if shortReadDepth and seqId in shortReadDepth:
                    meanDepth, normalizedDepth = shortReadDepth[seqId]
                    contigId += " %s %.01f %s %.2f"%(shortCovTags[0], meanDepth, shortCovTags[1], normalizedDepth)
                    contigCoverage = meanDepth
                if longReadDepth and seqId in longReadDepth:
                    meanDepth, normalizedDepth = longReadDepth[seqId]
                    contigId += " %s %.01f %s %.2f"%(longCovTags[0], meanDepth, longCovTags[1], normalizedDepth)
                    contigCoverage = max(meanDepth, contigCoverage)
                if "contigCircular" in details and contigIndex in details["contigCircular"]:
                    contigId += " circular=true"
                if len(seq) >= args.min_contig_length and contigCoverage >= args.min_contig_coverage:
                    OUT.write(contigId+"\n")
                    for i in range(0, len(seq), 60):
                        OUT.write(seq[i:i+60]+"\n")
                else:
                    suboptimalReads += contigId+"\n"+seq+"\n"
    if suboptimalReads:
        suboptimalReadsFile = re.sub("\..*", "_short_or_low_coverage.fasta", outputFile)
        with open(suboptimalReadsFile, "w") as SUBOPT:
            SUBOPT.write(suboptimalReads)
    if os.path.getsize(outputFile) < 10:
        LOG.write("failed to generate outputFile, return None\n")
        return None
    comment = "trimContigsByMinLength, input %s, output %s"%(inputFile, outputFile)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return outputFile

def runBandage(gfaFile, details):
    imageFormat = ".jpg"
    retval = None
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write("Bandage command =\n"+" ".join(command)+"\n")
        try:
            return_code = subprocess.call(command, shell=False, stderr=LOG)
            LOG.write("return code = %d\n"%return_code)
            if return_code == 0:
                retval = plotFile
            else:
                LOG.write("Error creating Bandage plot\n")
        except OSError as ose:
            comment = "Problem running Bandage: "+str(ose)
            LOG.write(comment+"\n")
            details['problem'].append(comment)
    return retval

def runUnicycler(details, threads=1, min_contig_length=0, prefix=""):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("runUnicycler: elapsed seconds = %f\n"%(time()-Start_time))
    command = ["unicycler", "-t", str(threads), "-o", '.']
    if min_contig_length:
        command.extend(("--min_fasta_length", str(min_contig_length)))
    command.extend(("--keep", "0")) # keep only assembly.gfa, assembly.fasta and unicycler.log
    command.append("--no_pilon")

    # put all read files on command line, let Unicycler figure out which type each is
    # apparently unicycler can only accept one read set in each class (I tried multiple ways to submit 2 paired-end sets, failed)
    short1 = None
    short2 = None
    unpaired = None
    long_reads = None
    for item in details['reads']:
        if 'superceded_by' in details['reads'][item]:
            continue
        files = details['reads'][item]['file'] # array of tuples (path, filename)
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
        command.extend(("--unpaired", " ".join(unpaired)))
    if long_reads:
        comand.extend(("--long", long_reads))
    # it is not quite right to send iontorrent data to spades through unicycler because the --iontorrent flag to spades will not be set

    LOG.write("Unicycler command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    unicyclerStartTime = time()
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

    details["assembly"] = { 
                'assembly_elapsed_time' : elapsedHumanReadable,
                'assembler': 'unicycler',
                'command_line': " ".join(command)
                }

    LOG.write("Duration of Unicycler run was %s\n"%(elapsedHumanReadable))

    if os.path.exists("unicycler.log"):
        unicyclerLogFile = prefix+"unicycler.log"
        shutil.move("unicycler.log", os.path.join(SaveDir, unicyclerLogFile))

    if not os.path.exists("assembly.fasta"):
        LOG.write("unicycler failed to generate assembly file.\n")
        details["problem"].append("unicycler failed to generate contigs file")
        return None
    details["contigCircular"] = []
    with open("assembly.fasta") as F:
        contigIndex = 1
        for line in F:
            if line.startswith(">"):
                if "circular=true" in line:
                    details["contigCircular"].append(contigIndex) 
                contigIndex += 1

    assemblyGraphFile = prefix+"assembly_graph.gfa"
    shutil.move("assembly.gfa", os.path.join(SaveDir, assemblyGraphFile))

    shutil.move("assembly.fasta", "contigs.fasta") #rename to canonical name
    return "contigs.fasta"

def runSpades(args, details):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("runSpades: elapsed seconds = %f\n"%(time()-Start_time))

    if args.illumina and args.iontorrent:
        comment = "SPAdes is not meant to process both Illumina and IonTorrent reads in the same run"
        details["problem"].append(comment)
        LOG.write(comment+"\n")
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
    if args.memory:
        command.extend(["-m", str(args.memory)])
    if args.recipe == "meta-spades":
        command.append("--meta")
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
    spadesStartTime = time()

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

    details["assembly"] = { 
                'assembly_elapsed_time' : elapsedHumanReadable,
                'assembler': 'spades',
                'command_line': " ".join(command)
                }

    LOG.write("Duration of SPAdes run was %s\n"%(elapsedHumanReadable))
    if not os.path.exists("contigs.fasta"):
        LOG.write("spades failed to generate contigs.fasta.\n")
        details["problem"].append("spades failed to generate contigs.fasta")
        return None
    try:
        spadesLogFile = args.prefix+"spades.log"
        shutil.move("spades.log", os.path.join(SaveDir, spadesLogFile))
        assemblyGraphFile = args.prefix+"assembly_graph.gfa"
        shutil.move("assembly_graph_with_scaffolds.gfa", os.path.join(SaveDir, assemblyGraphFile))
    except Exception as e:
        LOG.write(str(e))
    contigsFile = "contigs.fasta"
    return contigsFile

def runMinimap(contigFile, longReadFastq, details, threads=1, outformat='sam'):
    #LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    """
    Map long reads to contigs by minimap2 (read paf-file, readsFile; generate paf file).
    """
    LOG.write("runMinimap(%s, %s, details, %d, %s)\n"%(contigFile, longReadFastq, threads, outformat))
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

    #otherwise format as bam (and index)
    tempTime = time()
    sortThreads = max(int(threads/2), 1)
    fastqBase = longReadFastq
    fastqBase = re.sub(r"\..*", "", fastqBase)
    bamFile = contigFile+"_"+fastqBase+".bam"
    command = "samtools view -bS -@ %d %s | samtools sort -@ %d - -o %s"%(sortThreads, contigSam, sortThreads, bamFile)
    LOG.write("executing:\n"+command+"\n")
    return_code = subprocess.call(command, shell=True, stderr=LOG)
    LOG.write("samtools view|sort return code = %d, time=%d\n"%(return_code, time()-tempTime))
    os.remove(contigSam) #save a little space
    if return_code != 0:
        return None

    command = ["samtools", "index", bamFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools index return code = %d\n"%return_code)
    LOG.write('runMinimap returning %s\n'%bamFile)
    return bamFile

def runRacon(contigFile, longReadsFastq, details, threads=1):
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    Run racon on reads, read-to-contig-sam, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write('runRacon(%s, %s, details, %d)\n'%(contigFile, longReadsFastq, threads))
    readsToContigsSam = runMinimap(contigFile, longReadsFastq, details, threads, outformat='sam')
    raconStartTime = time()
    raconContigs = contigFile.replace(".fasta", ".racon.fasta")
    raconOut = open(raconContigs, 'w')
    command = ["racon", "-t", str(threads), "-u", longReadsFastq, readsToContigsSam, contigFile]
    tempTime = time()
    LOG.write("racon command: \n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL, stdout=raconOut)
    LOG.write("racon return code = %d, time = %d seconds\n"%(return_code, time()-raconStartTime))
    if return_code != 0:
        return None
    raconContigSize = os.path.getsize(raconContigs)
    if raconContigSize < 10:
        return None
    comment = "racon, input %s, output %s"%(contigFile, raconContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return raconContigs

def runBowtie(contigFile, shortReadFastq, details, threads=1, outformat='bam'):
    """
    index contigsFile, then run bowtie2, then convert sam file to pos-sorted bam and index
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    command = ["bowtie2-build", "--threads", str(threads),  contigFile, contigFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout and stderr to dev/null
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("bowtie2-build return code = %d\n"%return_code)
    if return_code != 0:
        return None

    command = ["bowtie2", "-p", str(threads), "-x", contigFile]
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

    sortThreads = max(int(threads/2), 1)
    bamFile = re.sub(r".sam$", ".bam", samFile)
    command = "samtools view -bS -@ %d %s | samtools sort -@ %d - -o %s"%(sortThreads, samFile, sortThreads, bamFile)
    LOG.write("executing:\n"+command+"\n")
    return_code = subprocess.call(command, shell=True, stderr=LOG)
    LOG.write("samtools return code = %d\n"%return_code)
    os.remove(samFile) #save a little space
    if return_code != 0:
        return None

    command = ["samtools", "index", bamFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools return code = %d\n"%return_code)
    return bamFile

def runPilon(contigFile, shortReadFastq, details, pilon_jar, threads=1):
    """ 
    polish contigs with short reads (illumina or iontorrent)
    first map reads to contigs with bowtie
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    if not os.path.exists(pilon_jar):
        comment = "jarfile %s not found when processing %s, giving up"%(pilon_jar, shortReadFastq)
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return

    bamFile = runBowtie(contigFile, shortReadFastq, details, threads=threads, outformat='bam')
    if not bamFile:
        return None
    command = ['java', '-Xmx32G', '-jar', pilon_jar, '--genome', contigFile]
    if ':' in shortReadFastq:
        command.extend(('--frags', bamFile))
    else:
        command.extend(('--unpaired', bamFile))
    pilonPrefix = contigFile.replace(".fasta", ".pilon")
    command.extend(('--outdir', '.', '--output', pilonPrefix, '--changes'))
    command.extend(('--threads', str(threads)))
    tempTime = time()
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("pilon return code = %d\n"%return_code)
    LOG.write("pilon duration = %d\n"%(time() - tempTime))
    if return_code != 0:
        return None
    pilonContigs = pilonPrefix+".fasta"
    numChanges = 0
    with open(pilonContigs.replace(".fasta", ".changes")) as CHANGES:
        for line in CHANGES:
            numChanges += 1
    LOG.write("Number of changes made by pilon was %d\n"%numChanges)

    comment = "pilon, input %s, output %s, num_changes = %d"%(contigFile, pilonContigs, numChanges)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return pilonContigs 

def calcReadDepth(bamfiles):
    """ Return dict of contig_ids to tuple of (coverage, normalized_coverage) """
    LOG.write("calcReadDepth(%s)\n"%" ".join(bamfiles))
    readDepth = {}
    command = ["samtools", "depth", "-a"]
    if type(bamfiles) is str:
        command.append(bamfiles)
    else:
        command.extend(bamfiles)
    LOG.write("command = "+" ".join(command)+"\n")
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    depthData = proc.communicate()[0]
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

    # calculate mean depth of contigs within "normal" boundary around overall mean
    lowerBound = totalMeanDepth * 0.5
    upperBound = totalMeanDepth * 1.5
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

def runCanu(details, threads=1, genome_size="5m", memory=250, prefix=""):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    canuStartTime = time()
    LOG.write("runCanu: elapsed seconds = %d\n"%(canuStartTime-Start_time))
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
    command = ["canu", "-d", '.', "-p", "canu", "useGrid=false", "genomeSize=%s"%genome_size]
    command.extend(["maxMemory=" + str(memory), "maxThreads=" + str(threads)])
    command.append("stopOnReadQuality=false")
    """
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    numLongReadFiles = 0
    if len(details['platform']['pacbio']):
        command.append("-pacbio-raw")
        command.extend(details['platform']['pacbio']) #allow multiple files
        numLongReadFiles += len(details['platform']['pacbio'])
    if len(details['platform']['fasta']):
        command.extend(details['platform']['fasta']) #allow multiple files
        numLongReadFiles += len(details['platform']['fasta'])
        comment = 'submitting fasta reads to canu, but calling them "pacbio": '+' '.join(details['platform']['fasta'])
        LOG.write(comment+"\n")
        details['problem'].append(comment)
    if len(details['platform']['nanopore']):
        command.append("-nanopore-raw")
        command.extend(details['platform']['nanopore']) #allow multiple files
        numLongReadFiles += len(details['platform']['nanopore'])
    if not numLongReadFiles:
        LOG.write("no long read files available for canu.\n")
        details["problem"].append("no long read files available for canu")
        return None
    LOG.write("canu command =\n"+" ".join(command)+"\n")

    canuStartTime = time()
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
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

    details["assembly"] = { 
                'assembly_elapsed_time' : elapsedHumanReadable,
                'assembler': 'canu',
                'command_line': " ".join(command)
                }

    LOG.write("Duration of canu run was %s\n"%(elapsedHumanReadable))
    if os.path.exists("canu.report"):
        shutil.move("canu.report", os.path.join(SaveDir, prefix+"canu_report.txt"))
    
    if not os.path.exists("canu.contigs.fasta"):
        LOG.write("canu failed to generate contigs file.\n")
        details["problem"].append("canu failed to generate contigs file")
        return None
    # rename to canonical contigs.fasta
    shutil.move("canu.contigs.fasta", "contigs.fasta")
    shutil.move("canu.contigs.gfa", os.path.join(SaveDir, prefix+"assembly_graph.gfa"))
    return "contigs.fasta"

def write_html_report(htmlFile, details):
    LOG.write("writing html report to %s\n"%htmlFile)
    HTML = open(htmlFile, 'w')
    HTML.write("<head><style>\n.a { text-indent: 50px }\n")
    HTML.write(".b {text-indent: 75px }\n")
    HTML.write("</style></head>")
    HTML.write("<h1>Genome Assembly Report</h1>\n")

    HTML.write("<h3>Input reads:</h3>\n")
    for item in details['reads']:
        if 'supercedes' in details['reads'][item]:
            continue # this is a derived item, not original input
        HTML.write(item+"<table class='a'>")
        for key in sorted(details['reads'][item]):
            if key == 'problems':
                continue
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
        HTML.write("</table>\n")
        if "problems" in details['reads'][item] and len(details['reads'][item]['problems']):
            HTML.write("<div class='b'><b>Issues with read set "+item+"</b>\n<ul>")
            for prob in details['reads'][item]['problems']:
                HTML.write("<li>"+prob+"\n")
            HTML.write("</ul></div>\n")
    
    if len(details["pre-assembly transformation"]):
        HTML.write("<h3>Pre-Assembly Transformations</h3>\n<div class='a'>\n")
        for line in details["pre-assembly transformation"]:
            HTML.write("<p>"+line+"\n")
        HTML.write("</div>\n")

    if "trim report" in details:
        HTML.write("<h3>Trimming Report</h3>\n<div class='a'>\n")
        for reads in details["trim report"]:
            HTML.write("<b>"+reads+"</b><ul>")
            for report in details["trim report"][reads]:
                if os.path.exists(os.path.join(SaveDir, report)):
                    HTML.write("<pre>\n")
                    HTML.write(open(os.path.join(SaveDir, report)).read())
                    HTML.write("\n</pre>\n")
        HTML.write("</div>\n")

    if len(details['derived_reads']):
        HTML.write("<h3>Transformed reads:</h3>\n")
        for item in details['derived_reads']:
            HTML.write(item+"<table class='a'>")
            for key in sorted(details['reads'][item]):
                if key in ('problems'):
                    continue
                HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
            HTML.write("</table>\n")
            if "problems" in details['reads'][item] and len(details['reads'][item]['problems']):
                HTML.write("<div class='b'><b>Issues with read set "+item+"</b>\n<ul>")
                for prob in details['reads'][item]['problems']:
                    HTML.write("<li>"+prob+"\n")
                HTML.write("</ul></div>\n")

    HTML.write("<h3>Assembly</h3>\n")
    if 'assembly' in details:
        HTML.write("<table class='a'>")
        for key in sorted(details['assembly']):
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['assembly'][key])))
        HTML.write("</table>\n")
    else:
        HTML.write("<p>None</p>\n")

    if "quast_txt" in details:
        HTML.write("<h3>Quast Report</h3>\n")
        HTML.write("<table class='a'>")
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_txt"], "Quast text report"))
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_html"], "Quast html report"))
        HTML.write("</table>\n")
        if os.path.exists(os.path.join(SaveDir, details["quast_txt"])):
            HTML.write("<pre>\n")
            HTML.write(open(os.path.join(SaveDir, details["quast_txt"])).read())
            HTML.write("\n</pre>\n")
    
    if len(details["post-assembly transformation"]):
        HTML.write("<h3>Post-Assembly Transformations</h3>\n<div class='a'>\n")
        for line in details["post-assembly transformation"]:
            HTML.write("<p>"+line+"<br>\n")
        HTML.write("</div>\n")

    if "Bandage plot" in details:
        path, imageFile = os.path.split(details["Bandage plot"])
        HTML.write("<h3>Bandage Plot</h3>\n")
        HTML.write("<div class='a'>")
        HTML.write("<img src='%s'>\n"%imageFile)
        HTML.write("</div>\n")
    HTML.close()


def main():
    global Start_time
    Start_time = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='p3_assembly')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    parser.add_argument('--sra', metavar='files', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('--recipe', choices=['unicycler', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'single-cell', 'auto'], help='assembler to use', default='auto')

    parser.add_argument('--racon_iterations', type=int, default=2, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=2, help='number of times to run pilon per short-read file', required=False)
    #parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--prefix', default='', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', default=Default_genome_size, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', default=300, help='save contigs of this length or longer', required=False)
    parser.add_argument('--min_contig_coverage', default=5, help='save contigs of this coverage or deeper', required=False)
    parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=250, help='RAM limit in Gb')
    parser.add_argument('--trim', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--pilon_jar', help='path to pilon executable or jar')
    parser.add_argument('--bandage', action='store_true', help='generate image of assembly path using Bandage')
    parser.add_argument('--params_json', help='JSON file with additional information.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    if args.params_json:
        parseJsonParameters(args)
    baseName = args.outputDirectory #"p3_assembly" 
    if len(args.prefix) > 0 and not args.prefix.endswith("_"):
        args.prefix += "_"
    global WorkDir
    WorkDir = baseName+"_work"
    if os.path.exists(WorkDir):
        shutil.rmtree(WorkDir)
    os.mkdir(WorkDir)
    WorkDir = os.path.abspath(WorkDir)
    global SaveDir
    SaveDir = os.path.abspath(os.path.join(WorkDir, "save"))
    if os.path.exists(SaveDir):
        shutil.rmtree(SaveDir)
    os.mkdir(SaveDir)
    logfileName = os.path.join(SaveDir, args.prefix + "p3_assembly.log")

    global LOG 
    sys.stderr.write("logging to "+logfileName+"\n")
    LOG = open(logfileName, 'w', 0) #unbuffered 
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(Start_time))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Work directory is "+WorkDir+"\n\n")
    LOG.write("Final output will be saved to "+SaveDir+"\n\n")
    details = { 'logfile' : logfileName }
    details["pre-assembly transformation"] = []
    details["post-assembly transformation"] = []
    details["original_items"] = []
    details["reads"] = {}
    details["problem"] = []
    details["derived_reads"] = []
    details["platform"] = {'illumina':[], 'iontorrent':[], 'pacbio':[], 'nanopore':[], 'fasta':[], 'anonymous':[]}

    if args.illumina:
        platform='illumina'
        for item in args.illumina:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved)

    if args.iontorrent:
        platform='iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved)

    if args.pacbio:
        for item in args.pacbio:
            registerReads(item, details, platform='pacbio')

    if args.nanopore:
        for item in args.nanopore:
            registerReads(item, details, platform='nanopore')

    if args.sra:
        fetch_sra_files(args.sra, details)

    if args.anonymous_reads:
        categorize_anonymous_read_files(args, details)

    # move into working directory so that all files are local
    original_working_directory = os.getcwd()
    os.chdir(WorkDir)

    if args.trim and len(details['platform']['illumina'] + details['platform']['iontorrent']):
        trimGalore(details, threads=args.threads)

    if args.recipe == "auto":
        #now must decide which assembler to use
        if True:
            # original rule: if any illumina or iontorrent reads present, use Unicycler (long-reads can be present), else use canu for long-reads
            if len(details['platform']['illumina'] + details['platform']['iontorrent']):
                args.recipe = "unicycler"
            else:
                args.recipe = "canu"
        else:
            # alternative rule: if any long reads present, use canu
            if len(details['platform']['pacbio'] + details['platform']['nanopore']):
                args.recipe = "canu"
            else:
                args.recipe = "unicycler"
    if "spades" in args.recipe or args.recipe == "single-cell":
        contigs = runSpades(args, details)
    elif args.recipe == "unicycler":
        contigs = runUnicycler(details, threads=args.threads, min_contig_length=args.min_contig_length, prefix=args.prefix)
    elif args.recipe == "canu":
        contigs = runCanu(details, threads=args.threads, genome_size=args.genome_size, memory=args.memory, prefix=args.prefix)
    else:
        LOG.write("cannot interpret args.recipe: "+args.recipe)

    if contigs and os.path.getsize(contigs):
        # now run racon with each long-read file
        for i in range(0, args.racon_iterations):
            for longReadFile in details['reads']:
                if details['reads'][longReadFile]['length_class'] == 'long':
                    LOG.write("runRacon(%s, %s, details, threads=%d)\n"%(contigs, longReadFile, args.threads))
                    raconContigFile = runRacon(contigs, longReadFile, details, threads=args.threads)
                    if raconContigFile is not None:
                        contigs = raconContigFile
                    else:
                        break
        
    if contigs and os.path.getsize(contigs):
        # now run pilon with each short-read file
        for iteration in range(0, args.pilon_iterations):
            numChanges = 0
            for shortReadFastq in details['reads']:
                if 'superceded_by' in details['reads'][shortReadFastq]:
                    continue
                if details['reads'][shortReadFastq]['length_class'] == 'short':
                    LOG.write("runPilon(%s, %s, details, %s, threads=%d) iteration=%d\n"%(contigs, shortReadFastq, args.pilon_jar, args.threads, iteration))
                    pilonContigFile = runPilon(contigs, shortReadFastq, details, args.pilon_jar, threads=args.threads)
                    if pilonContigFile is not None:
                        contigs = pilonContigFile
                    else:
                        break
            if not numChanges:
                break
        
    if contigs and os.path.getsize(contigs):
        shortReadDepth={}
        longReadDepth={}
        bamFiles=[]
        for reads in details['reads']:
            if details['reads'][reads]['length_class'] == 'short':
                bam = runBowtie(contigs, reads, details, threads=args.threads, outformat='bam')
                if bam:
                    bamFiles.append(bam)
        if len(bamFiles):
            shortReadDepth = calcReadDepth(bamFiles)
        bamFiles=[]
        for reads in details['reads']:
            if details['reads'][reads]['length_class'] == 'long':
                bam = runMinimap(contigs, reads, details, threads=args.threads, outformat='bam')
                if bam:
                    bamFiles.append(bam)
        if len(bamFiles):
            longReadDepth = calcReadDepth(bamFiles)
        saveContigsFile = os.path.join(SaveDir, args.prefix+"contigs.fasta")
        filteredContigs = filterContigsByMinLength(contigs, args, details, shortReadDepth=shortReadDepth, longReadDepth=longReadDepth)
        if filteredContigs:
            contigs = filteredContigs
    if contigs and os.path.getsize(contigs):
        runQuast(contigs, args, details)
        shutil.move(contigs, os.path.join(SaveDir, args.prefix+"contigs.fasta"))

    gfaFile = os.path.join(SaveDir, args.prefix+"assembly_graph.gfa")
    if os.path.exists(gfaFile):
        bandagePlot = runBandage(gfaFile, details)
        details["Bandage plot"] = bandagePlot

    htmlFile = os.path.join(SaveDir, args.prefix+"assembly_report.html")
    write_html_report(htmlFile, details)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - Start_time)/3600))
    LOG.close()
    fp = file(os.path.join(SaveDir, args.prefix+"run_details.txt"), "w")
    json.dump(details, fp, indent=2, sort_keys=True)
    fp.close()


if __name__ == "__main__":
    main()
