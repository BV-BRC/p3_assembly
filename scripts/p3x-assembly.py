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

def registerReads(details, reads, platform=None, interleaved=False):
    """
    create an entry in details for these reads
    move read files to working directory to allow relative paths
    """
    if "reads" not in details:
        details["reads"] = {}
    if 'original_items' not in details:
        details['original_items'] = []
    if reads in details["original_items"]:
        comment = "dulicate registration of reads %s"%reads
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    details['original_items'].append(reads)
    
    readStruct = {}
    readStruct["file"] = []
    readStruct["path"] = []
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
        file1 = os.path.split(read1)[1]
        file2 = os.path.split(read2)[1]
        os.symlink(os.path.abspath(read1), os.path.join(WorkDir, file1))
        os.symlink(os.path.abspath(read2), os.path.join(WorkDir, file2))
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
        file1 = os.path.split(reads)[1]
        os.symlink(os.path.abspath(reads), os.path.join(WorkDir, file1))
        readStruct["file"].append(file1)
        readStruct["path"].append(reads)
        registeredName = file1

    details['reads'][registeredName]=readStruct
    if platform:
        readStruct['platform'] = platform
        if platform not in details["platform"]:
            details["platform"][platform] = []
        details["platform"][platform].append(registeredName)
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
                if pe["platform"] == 'illumina':
                    if not args.illumina:
                        args.illumina = []
                    args.illumina.append(":".join((data["read1"], data["read2"]))) 
                elif pe["platform"] == 'iontorrent':
                    if not args.iontorrent:
                        args.iontorrent = []
                    args.iontorrent.append(":".join((data["read1"], data["read2"]))) 
                else:
                    if not args.anonymous_reads:
                        args.anonymous_reads = []
                    args.anonymous_reads.append(":".join((data["read1"], data["read2"]))) 
        if "single_end_libs" in data:
            for pe in data["single_end_libs"]:
                if pe["platform"] == 'illumina':
                    if not args.illumina:
                        args.illumina = []
                    args.illumina.append(data["read"]) 
                elif pe["platform"] == 'iontorrent':
                    if not args.iontorrent:
                        args.iontorrent = []
                    args.iontorrent.append(data["read"]) 
                elif pe["platform"] == 'pacbio':
                    if not args.pacbio:
                        args.pacbio = []
                    args.pacbio.append(data["read"]) 
                elif pe["platform"] == 'nanopore':
                    if not args.nanopore:
                        args.nanopore = []
                    args.nanopore.append(data["read"]) 
                else:
                    if not args.anonymous_reads:
                        args.anonymous_reads = []
                    args.anonymous_reads.append(":".join((data["read1"], data["read2"]))) 
        if "srr_ids" in data:
            if not args.sra:
                args.sra=[]
                args.sra.extend(data["srr_ids"])

def studyPairedReads(item, details):
    func_start = time()
    LOG.write("studyPairedReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start - Start_time))
    if "reads" not in details:
        details["reads"] = {}
    if item not in details["reads"]:
        details["reads"][item] = {}
    details['reads'][item]['layout'] = 'paired-end'
    details['reads'][item]['avg_len'] = 0
    details['reads'][item]['num_reads'] = 0
    details['reads'][item]['problems'] = []
    file1, file2 = item.split(":")
    if file1.endswith("gz"):
        F1 = gzip.open(file1)
        F2 = gzip.open(file2)
    else:
        F1 = open(file1)
        F2 = open(file2)

    line = F1.readline()
    sample_read_id = line.split(' ')[0]
    F1.seek(0)
    if sample_read_id.startswith('>'):
        details['reads'][item]['platform']='fasta'
        if 'fasta' not in details['platform']:
            details['platform']['fasta'] = []
        details['platform']['fasta'].append(item)
        #studyFastaReads(F1)
        #studyFastaReads(F2)
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
            if seqQualLenMatch and (seqLen1 != len(line1)-1 or seqLen2 != len(line2)-1):
                seqQualLenMatch = False
                details['reads'][item]["problems"].append("sequence and quality strings differ in length at read %d"%readNumber)
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
    for line in F:
        if i % 4 == 0 and not sample_read_id:
            sample_read_id = line.split(' ')[0] # get part up to first space, if any 
        elif i % 4 == 1:
            seqLen = len(line)-1
        elif i % 4 == 3:
            if seqQualLenMatch and (seqLen != len(line)-1):
                seqQualLenMatch = False
                details['reads'][item]["problems"].append("sequence and quality strings differ in length at read %d"%readNumber)
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

    LOG.write("duration of studySingleReads was %d seconds\n"%(time() - func_start))
    return

def studyFastaReads(F, details):
    """
    assume format is fasta
    F is opened file
    count reads, calc total length, mean, max, min
    """
    seq = ""
    seqLen = 0
    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    i = 0
    for line in F:
        if line.startswith(">"):
            if seq:
                seqLen = len(seq)
                totalReadLength += seqLen
                maxReadLength = max(maxReadLength, seqLen) 
                minReadLength = min(minReadLength, seqLen)
                minQualScore = min(minQualScore + line.rstrip())
                maxQualScore = max(maxQualScore + line.rstrip())
                readNumber += 1
                seq = ""
        else:
            seq += line.rstrip()
    if seq:
        seqLen = len(seq)
        totalReadLength += seqLen
        maxReadLength = max(maxReadLength, seqLen) 
        minReadLength = min(minReadLength, seqLen)
        minQualScore = min(minQualScore + line.rstrip())
        maxQualScore = max(maxQualScore + line.rstrip())
        readNumber += 1
        seq = ""

    avgReadLength = totalReadLength/readNumber
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 
    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, avgReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][avgReadLength >= Max_short_read_length]

    LOG.write("duration of studySingleReads was %d seconds\n"%(time() - func_start))
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
            return "sra" # 
# NOTE: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
    if re.match(r"@\S+/\S+/\S+_\S+$", read_id): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
        return "pacbio_ccs_subread" # 
    if re.match(r"@\S+/\S+$", read_id): #@<MovieName>/<ZMW_number> 
        return "pacbio_ccs" # 
#@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
    if re.match(r"@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", read_id): #based on one example, need to test more 
        return "nanopore" # 
    return "na"

def trimGalore(args, details):
    startTrimTime = time()
    LOG.write("\ntrimGalore() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    if "trim report" not in details:
        details["trim report"] = {}
    for platformReads in (args.illumina, args.iontorrent):
        for index, reads in enumerate(platformReads): 
            command = ['trim_galore', '-j', str(args.threads), '-o', '.']
            if ':' in reads:
                splitReads = reads.split(":")
                #-j 4 -o testTrim --length 30 --paired SRR1395326_1_10pct.fastq SRR1395326_2_10pct.fastq
                command.extend(["--paired", splitReads[0], splitReads[1]])
                LOG.write("command: "+" ".join(command))
                return_code = subprocess.call(command, shell=False, stderr=LOG)
                LOG.write("return code = %d\n"%return_code)
                trimReads1 = os.path.basename(splitReads[0])
                trimReport1 = trimReads1+"_trimming_report.txt"
                trimReads1 = trimReads1.replace(".fastq", "_val_1.fq")
                trimReads2 = os.path.basename(splitReads[1])
                trimReport2 = trimReads2+"_trimming_report.txt"
                trimReads2 = trimReads2.replace(".fastq", "_val_2.fq")
                comment = "trim_galore, input %s, output %s"%(reads, trimReads1+":"+trimReads2)
                if os.path.exists(trimReads1) and os.path.exists(trimReads2):
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                    platformReads[index] = trimReads1+":"+trimReads2
                    if os.path.exists(trimReport1) and os.path.exists(trimReport2):
                        shutil.move(trimReport1, os.path.join(SaveDir, trimReport1))
                        shutil.move(trimReport2, os.path.join(SaveDir, trimReport2))
                        details["trim report"][reads]=[trimReport1, trimReport2]
                else:
                    comment = "Problem during trim_galore: expected files not found: "+trimReads1 +", "+trimReads2
                    details["problem"].append(comment)
                    LOG.write(comment)
            else:
                command.append(reads)
                LOG.write("command: "+" ".join(command))
                return_code = subprocess.call(command, shell=False, stderr=LOG)
                LOG.write("return code = %d\n"%return_code)
                trimReads = os.path.basename(reads)
                trimReport = trimReads+"_trimming_report.txt"
                trimReads = trimReads.replace(".fastq", "_trimmed.fq")
                comment = "trim_galore, input %s, output %s"%(reads, trimReads)
                if os.path.exists(trimReads):
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                    platformReads[index] = trimReads
                    if os.path.exists(trimReport):
                        shutil.move(trimReport, os.path.join(SaveDir, trimReport))
                        details["trim report"][reads]=[trimReport]
                else:
                    comment = "Problem during trim_galore: expected files not found: "+trimReads
                    LOG.write(comment)
                    details["problem"].append(comment)
    LOG.write("trim_galore trimming completed, duration = %d seconds\n\n\n"%(time()-startTrimTime))

def determineReadType(filename, details=None):
    srf_time = time()
    LOG.write("determineReadType() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(srf_time)), srf_time-Start_time))
    # figures out Read_file_type
    #return read_file_type and sample of read ids
    read_file_type = 'na'
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
        comment = "in determineReadType for %s: text sample (length %d) lacks at least 2 lines"%(filename, len(text))
        LOG.write(comment+"\n")
        details.problem.append(comment)
    if lines[0].startswith("@"):
        read_file_type = 'fastq'
    elif lines[0].startswith(">"):
        read_file_type = 'fasta'
    if read_format == 'fastq':
        avg_read_length = 0
        readLengths = []
        for i, line in enumerate(lines):
            if i % 4 == 0:
                read_id_sample.append(line.split(' ')[0]) # get part up to first space, if any 
            elif i % 4 == 1:
                readLengths.append(len(line)-1)
        avg_read_length = sum(readLengths)/float(len(readLengths))
        read_file_type = inferPlatform(read_id_sample[0], avg_read_length)
        LOG.write("fastq read type %s average read length %.1f\n"%(read_file_type, avg_read_lengt))
    if details:
        if not "inferred_read_type" in details:
            details["inferred_read_type"] = {}
        details["inferred_read_type"][filename] = Read_file_type[filename]
        if not "estimated_read_length" in details:
            details["estimated_read_length"] = {}
        details["estimated_read_length"][filename] = Avg_read_length[filename]
    LOG.write("duration of determineReadType was %d seconds\n"%(time() - srf_time))
    return read_file_type, read_id_sample

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
            read_file_type[item], read_id_sample[item] = determineReadType(item, details)
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
            read_file_type[filename1], read_id_sample[filename1] = determineReadType(filename1, details)
            comment = "interpreting %s type as %s"%(filename1, read_file_type[filename1])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)
        if filename2 not in read_file_type:
            read_file_type[filename2], read_id_sample[filename2] = determineReadType(filename2, details)
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
            valid_pairs.append(item)
            membersOfPairs.add(filename1)
            membersOfPairs.add(filename2)
        else: #move over to single files
            valid_singles.add(filename1)
            valid_singles.add(filename2)
    
    # some items on singleFiles may not be valid (may be paired up)
    for item in singleFiles:
        if item not in membersOfPairs:
            valid_singles.add(item)

    for item in valid_pairs + valid_singles:
        registerReads(details, item, platform=read_file_type[item], interleaved = args.interleaved and item in args.interleaved)

    return

def get_runinfo(run_accession, log=None):
    """ take sra run accession (like SRR123456)
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    Altered from versionin sra_tools to handle case of multiple sra runs returned by query.
    If efetch doesn't work, try scraping the web page.
    """
    LOG.write("get_runinfo(%s)\n"%run_accession)
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
        if not runinfo:
            return None
        if runinfo['Platform'] not in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
            # rescue case of mis-alignment between keys and values
            if log:
                log.write("problem in get_runinfo: sra.cgi returned:\n"+text+"\n")
            for val in values:
                if val in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
                    runinfo['Platform'] = val
                    break
        if runinfo['LibraryLayout'] not in ('PAIRED', 'SINGLE'):        
            if log:
                log.write("Need to search for LibraryLayout: bad value: %s\n"%runinfo['LibraryLayout'])
            for val in values:
                if val in ('PAIRED', 'SINGLE'):
                    runinfo['LibraryLayout'] = val
                    break

    if not runinfo:
        if log:
            log.write("Problem, normal runinfo request failed. Trying alternative from web page.\n")
        runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run="+run_accession
        text = urllib2.urlopen(runinfo_url).read()
        if re.search("<td>Illumina", text, re.IGNORECASE):
            runinfo = {'Platform': 'ILLUMINA'}
        elif re.search("<td>PacBio", text, re.IGNORECASE):
            runinfo = {'Platform': 'PACBIO_SMRT'}
        elif re.search("<td>Oxford", text, re.IGNORECASE):
            runinfo = {'Platform': 'OXFORD_NANOPORE'}
        elif re.search("<td>Ion Torrent", text, re.IGNORECASE):
            runinfo = {'Platform': 'ION_TORRENT'}
    return runinfo

def fetch_one_sra(sra, run_info=None, log=sys.stderr):
    """ requires run_info to know which program to use
    """
    if not run_info:
        run_info = get_runinfo(sra, log)
    programToUse = "fasterq-dump" # but not appropriate for pacbio or nanopore
    if runinfo['Platform'].startswith("PACBIO") or runinfo['Platform'].startsWith("OXFORD_NANOPORE"):
        programToUse = 'fastq-dump'
    command = [programToUse, "--split-files", sra]
    log.write("command = "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=log)
    LOG.write("return_code = %d, time=%d seconds\n"%(return_code, time()-stime))
    if return_code != 0:
        log.write("Problem, %s return code was %d\n"%(programToUse, return_code))

        log.write("Try one more time.\n")
        return_code = subprocess.call(command, shell=False, stderr=LOG)
        log.write("Return code on second try was %d\n"%return_code)
        if return_code != 0:
            LOG.write("Giving up on %s\n"%sra)
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
        runinfo = get_runinfo(sra)
        if not runinfo:
            LOG.write("runinfo for %s was empty, giving up\n"%sra)
            continue
        LOG.write("Runinfo for %s reports platform = %s and LibraryLayout = %s\n"%(sra, runinfo["Platform"], runinfo['LibraryLayout']))

        fetch_one_sra(sraFull, runinfo, LOG)
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
                comment = "for library %s, number of files was %s"%(sra, len(fastqFiles))
                details['problem'].append(comment)
                LOG.write(comment+"\n")
                subprocess.run("cat %s*fastq > %s.fastq"%(sra, sra), shell=True)
                item = sra+".fastq"
            if not item:
                comment = "for %s no fastq file found"%sra
                details['problem'].append(comment)
                LOG.write(comment+"\n")
            continue # failed on that sra
    
        if runinfo["Platform"] == "ILLUMINA":
            platform = "illumina"
        elif runinfo["Platform"] == "IONTORRENT":
            platform = "iontorrent"
        elif runinfo["Platform"] == "PACBIO":
            platform = "pacbio"
        elif runinfo["Platform"] == "OXFORD_NANOPORE":
            platform = "nanopore"
        else:
            platform="anonymous"
        registerReads(details, item, platform=platform)

    return

def writeSpadesYamlFile(details):
    LOG.write("writeSpadesYamlFile: elapsed seconds = %f\n"%(time()-Start_time))
    outfileName = "spades_yaml_file.txt"
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    for platorm in details['platform']:
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
            paired_end_reads[0].append(f)
        elif "%" in item:
            f = details['reads'][item]['file'][0]
            mate_pair_reads[0].append(f)
            f = details['reads'][item]['file'][1]
            mate_pair_reads[0].append(f)
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
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(paired_end_reads[0]))
        OUT.write("\"\n    ],\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(paired_end_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement=True
    if len(mate_pair_reads[0]):
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"rf\",\n")
        OUT.write("    type: \"mate-pairs\",\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[0]))
        OUT.write("\"\n    ]\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement=True
    if 'pacbio' in details['platform']:
        pacbio_reads = []
        for item in details['platform']['pacbio']:
            f = details['reads'][item]['file'][0]
            pacbio_reads.append(os.path.abspath(f))
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if 'nanopore' in details['platform']:
        nanopore_reads = []
        for item in details['platform']['nanopore']:
            f = details['reads'][item]['path'][0]
            nanopore_reads.append(os.path.abspath(f))
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if 'fasta' in details['platform']:
        fasta_reads = []
        for item in details['platform']['fasta']:
            f = details['reads'][item]['file'][0]
            nanopore_reads.append(os.path.abspath(f))
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"untrusted-contigs\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
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
    details["quast_txt"] = "quast_report.txt"
    details["quast_html"] = "quast_report.html"

def filterContigsByMinLength(inputFile, outputFile, args, details, shortReadDepth=None, longReadDepth=None):
    """ 
    Write only sequences at or above min_length to output file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("filterContigsByMinLength(%s, %s, %d) elapsed seconds = %f\n"%(inputFile, outputFile, args.min_contig_length, time()-Start_time))
    with open(inputFile) as IN:
        with open(outputFile, 'w') as OUT:
            shortCovTags = ['cov', 'norm_cov']
            #shortCovTags = [['cov', 'norm_cov'], ['cov', 'norm_cov']][longReadDepth != None]
            longCovTags = ['cov_long', 'norm_cov_long']
            seqId=None
            seq = ""
            i = 1
            for line in IN:
                m = re.match(">(\S+)", line)
                if m:
                    if len(seq) >= args.min_contig_length:
                        contigId = ">"+args.prefix+"contig_%d length %5d"%(i, len(seq))
                        if shortReadDepth and seqId in shortReadDepth:
                            meanDepth, length, normalizedDepth = shortReadDepth[seqId]
                            contigId += " %s %.01f %s %.2f"%(shortCovTags[0], meanDepth, shortCovTags[1], normalizedDepth)
                        if longReadDepth and seqId in longReadDepth:
                            meanDepth, length, normalizedDepth = longReadDepth[seqId]
                            contigId += " %s %.01f %s %.2f"%(longCovTags[0], meanDepth, longCovTags[1], normalizedDepth)
                        if "contigCircular" in details and seqId in details["contigCircular"]:
                            contigId += " circular=true"
                        OUT.write(contigId+"\n")
                        for i in range(0, len(seq), 60):
                            OUT.write(seq[i:i+60]+"\n")
                        i += 1
                    seq = ""
                    seqId = m.group(1)
                else:
                    seq += line.rstrip()
            if len(seq) >= args.min_contig_length:
                contigId = ">"+args.prefix+"contig_%d length %5d"%(i, len(seq))
                if shortReadDepth and seqId in shortReadDepth:
                    meanDepth, length, normalizedDepth = shortReadDepth[seqId]
                    contigId += " %s %.01f %s %.2f"%(shortCovTags[0], meanDepth, shortCovTags[1], normalizedDepth)
                    if int(float(length)) != len(seq):
                        LOG.write("len for %s conflict: %d vs %d\n"%(seqId, int(float(length)), len(seq)))
                if longReadDepth and seqId in longReadDepth:
                    meanDepth, length, normalizedDepth = longReadDepth[seqId]
                    contigId += " %s %.01f %s %.2f"%(longCovTags[0], meanDepth, longCovTags[1], normalizedDepth)
                OUT.write(contigId+"\n")
                for i in range(0, len(seq), 60):
                    OUT.write(seq[i:i+60]+"\n")
    comment = "trimContigsByMinLength, input %s, output %s"%(inputFile, outputFile)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)

def runBandage(args, details):
    gfaFile = os.path.join(SaveDir, args.prefix+"assembly_graph.gfa")
    imageFormat = ".jpg"
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write("Bandage command =\n"+" ".join(command)+"\n")
        try:
            return_code = subprocess.call(command, shell=False, stderr=LOG)
            LOG.write("return code = %d\n"%return_code)
            if return_code == 0:
                details["Bandage plot"] = args.prefix+"assembly_graph.plot"+imageFormat
            else:
                LOG.write("Error creating Bandage plot\n")
        except OSError as ose:
            comment = "Problem running Bandage: "+str(ose)
            LOG.write(comment+"\n")
            details['problem'].append(comment)


def runUnicycler(args, details):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("runUnicycler: elapsed seconds = %f\n"%(time()-Start_time))
    command = ["unicycler", "-t", str(args.threads), "-o", '.']
    #, "--pilon_path", args.pilon_jar]
    if args.min_contig_length:
        command.extend(("--min_fasta_length", str(args.min_contig_length)))
    command.extend(("--keep", "0")) # keep only assembly.gfa, assembly.fasta and unicycler.log
    #if args.no_pilon:
    command.append("--no_pilon")

    # put all read files on command line, let Unicycler figure out which type each is
    for item in details['reads']:
        files = details['reads'][item]['file'] # array of tuples (path, filename)
        if len(files) > 1:
            command.extend(("-1", files[0], "-2", files[1]))
        else:
            if details['reads'][item]['length_class'] == 'short':
                size_tag = '-s'
            else:
                size_tag = '-l'
            command.extend((size_tag, files[0]))
    # it is not quite right to send iontorrent data to spades through unicycler because the --iontorrent flag to spades will not be set

    LOG.write("Unicycler command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    unicyclerStartTime = time()

    return_code = subprocess.call(command, shell=False, stdout=LOG)
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
                'command_line': command,
                }

    LOG.write("Duration of Unicycler run was %s\n"%(elapsedHumanReadable))

    unicyclerLogFile = args.prefix+"unicycler.log"
    shutil.move("unicycler.log", os.path.join(SaveDir, unicyclerLogFile))

    if not os.path.exists("assembly.fasta"):
        LOG.write("unicycler failed to generate assembly file.\n")
        details["problem"].append("unicycler failed to generate contigs file")
        return None
    details["contigCircular"] = []
    with open("assembly.fasta") as F:
        for line in F:
            if line.startswith(">"):
                details["contigCircular"].append("circular=true" in line) # append True or False

    assemblyGraphFile = args.prefix+"assembly_graph.gfa"
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
    if args.singlecell:
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

    return_code = subprocess.call(command, shell=False, stderr=LOG)
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
                'command_line': command,
                }

    LOG.write("Duration of SPAdes run was %s\n"%(elapsedHumanReadable))
    if not os.path.exists("contigs.fasta"):
        LOG.write("spades failed to generate contigs.fasta")
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

def runMinimap(contigFile, longReadFastq, args, details):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    1. Index contigs for minimap2 (generate mmi file).
    2. Map long reads to contigs by minimap2 (read paf-file, readsFile; generate paf file).
    3. Run racon on reads, read-to-contig-map, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    # index contig sequences
    contigIndex = contigFile.replace(".fasta", ".mmi")
    command = ["minimap2", "-t", str(args.threads), "-d", contigIndex, contigFile] 
    tempTime = time() 
    LOG.write("minimap2 index command:\n"+' '.join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("minimap2 index return code = %d\n"%return_code)
    if return_code != 0:
        return None
    details['minimap2_index_time'] = time() - tempTime

    # map long reads to contigs
    contigSam = contigFile.replace(".fasta", ".sam")
    command = ["minimap2", "-t", str(args.threads), "-a", "-o", contigSam, contigFile, longReadFastq]
    tempTime = time()
    LOG.write("minimap2 map command:\n"+' '.join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("minimap2 map return code = %d\n"%return_code)
    if return_code != 0:
        return None
    details['minimap2_map_time'] = time() - tempTime
    return contigSam

def runRacon(contigFile, longReadsFastq, args, details):
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    Run racon on reads, read-to-contig-sam, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    readsToContigsSam = runMinimap(contigFile, longReadsFastq, args.threads, details)
    raconStartTime = time()
    raconContigs = contigFile.replace(".fasta", ".racon.fasta")
    raconOut = open(raconContigs, 'w')
    command = ["racon", "-t", str(args.threads), "-u", longReadsFastq, readsToContigsSam, contigFile]
    tempTime = time()
    LOG.write("racon command: \n"+' '.join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG, stdout=raconOut)
    LOG.write("racon return code = %d\n"%return_code)
    LOG.write("Racon total time = %d seconds\n"%(time()-raconStartTime))
    if return_code != 0:
        return None
    details['racon time'] = time()-tempTime
    raconContigSize = os.path.getsize(raconContigs)
    if raconContigSize < 10:
        return None
    comment = "racon, input %s, output %s"%(contigFile, raconContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return raconContigs

def runBowtie(contigFile, shortReadFastq, threads=1, output='bam'):
    """
    index contigsFile, then run bowtie2, then convert sam file to pos-sorted bam and index
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    command = ["bowtie2-build", "--threads", str(threads),  contigFile, contigFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
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
    fastqBase = re.sub("\..*", "", fastqBase)
    samFile = contigFile+"_"+fastqBase+".sam"
    command.extend(('-S', samFile))
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("bowtie2 return code = %d\n"%return_code)
    if return_code != 0:
        return None
    if output == 'sam':
        return samFile

    sortThreads = max(int(threads/2), 1)
    bamFile = re.sub(".sam$", ".bam", samFile)
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

def runPilon(contigFile, shortReadFastq, args, details):
    """ 
    polish contigs with short reads (illumina or iontorrent)
    first map reads to contigs with bowtie
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    if not args.pilon_jar:
        comment = "pilon_jar not defined when processing %s, giving up"%shortReadFastq
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return

    bamFile = runBowtie(contigFile, shortReadFastq, args.threads)
    if not bamFile:
        return None
    command = ['java', '-Xmx32G', '-jar', args.pilon_jar, '--genome', contigFile]
    if ':' in shortReadFastq:
        command.extend(('--frags', bamFile))
    else:
        command.extend(('--unpaired', bamFile))
    pilonPrefix = contigFile.replace(".fasta", ".pilon")
    command.extend(('--outdir', '.', '--output', pilonPrefix, '--changes'))
    command.extend(('--threads', str(args.threads)))
    tempTime = time()
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
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

    comment = "pilon, input %s, output %s"%(contigFile, pilonContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return pilonContigs 

def calcReadDepth(bamfiles):
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
                readDepth[prevContig] = [meanDepth, length, 0]
                depthSum = 0
                length = 0
            prevContig = contig
        totalLength += 1
        length += 1
        depthSum += depth
        totalDepthSum += depth
    if prevContig is not None:
        meanDepth = depthSum/length
        readDepth[prevContig] = [meanDepth, length, 0]

    LOG.write("total length for depth data = %d\n"%totalLength)
    LOG.write("total depth = %.1f\n"%totalDepthSum)
    meanDepth = 0
    if length > 0:
        meanDepth = depthSum/length
    readDepth[prevContig] = [meanDepth, length, 0]
    LOG.write("len(readDepth) = %d\n"%len(readDepth))
    totalMeanDepth = 0
    if totalLength > 0:
        totalMeanDepth = totalDepthSum/totalLength
    oneXSum = 0.0
    oneXLen = 0
    for c in readDepth:
        meanDepth, length = readDepth[c][:2]
        if meanDepth > totalMeanDepth * 0.5 and meanDepth < totalMeanDepth * 1.5:
            oneXSum += meanDepth * length
            oneXLen += length
    oneXDepth = 1
    if oneXLen > 0 and oneXSum > 0:
        oneXDepth = oneXSum/oneXLen
    for c in readDepth:
        meanDepth, length = readDepth[c][:2]
        normalizedDepth = meanDepth / oneXDepth
        readDepth[c][2] = normalizedDepth
    LOG.write("len(readDepth) = %d\n"%len(readDepth))
    return readDepth

def runMinimap(contigFile, longReadFastq, threads=1, output = 'sam'):
    """
    1. Index contigs for minimap2 (generate mmi file).
    2. Map long reads to contigs by minimap2, generate sam file.
    3. IF output specifies 'sam', then return name of sam file.
    4. Otherwise, convert to bam, index, and return name of bam file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    # index contig sequences
    contigIndex = contigFile.replace(".fasta", ".mmi")
    command = ["minimap2", "-t", str(threads), "-d", contigIndex, contigFile] 
    tempTime = time() 
    LOG.write("minimap2 index command:\n"+' '.join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("minimap2 index return code = %d\n"%return_code)
    if return_code != 0:
        return None

    # map long reads to contigs
    contigSam = contigFile.replace(".fasta", ".sam")
    command = ["minimap2", "-t", str(threads), "-a", "-o", contigSam, contigFile, longReadFastq]
    tempTime = time()
    LOG.write("minimap2 map command:\n"+' '.join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("minimap2 map return code = %d\n"%return_code)
    if return_code != 0:
        return None
    LOG.write("minimap2_map_time = %d\n"%(time()-tempTime))

    if output == 'sam':
        return contigSam

    #otherwise format as bam (and index)
    tempTime = time()
    sortThreads = max(int(threads/2), 1)
    fastqBase = longReadFastq
    fastqBase = re.sub("\..*", "", fastqBase)
    bamFile = contigFile+"_"+fastqBase+".bam"
    command = "samtools view -bS -@ %d %s | samtools sort -@ %d - -o %s"%(sortThreads, contigSam, sortThreads, bamFile)
    LOG.write("executing:\n"+command+"\n")
    return_code = subprocess.call(command, shell=True, stderr=LOG)
    LOG.write("samtools return code = %d, time=%d\n"%(return_code, time()-tempTime))
    os.remove(contigSam) #save a little space
    if return_code != 0:
        return None

    command = ["samtools", "index", bamFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools return code = %d\n"%return_code)
    return bamFile

def runCanu(args, details):
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
    command = ["canu", "-d", '.', "-p", "canu", "useGrid=false", "genomeSize=%s"%args.genome_size]
    command.extend(["maxMemory=" + str(args.memory), "maxThreads=" + str(args.threads)])
    command.append("stopOnReadQuality=false")
    """
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    if 'pacbio' in details['platform']:
        command.append("-pacbio-raw")
        command.extend(details['platform']['pacbio']) #allow multiple files
    if 'nanopore' in details['platform']:
        command.append("-nanopore-raw")
        command.extend(details['platform']['nanopore']) #allow multiple files
    LOG.write("canu command =\n"+" ".join(command)+"\n")

    canuStartTime = time()
    return_code = subprocess.call(command, shell=False, stdout=LOG)
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
                'command_line': command,
                }

    LOG.write("Duration of canu run was %s\n"%(elapsedHumanReadable))
    
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
        HTML.write(item+"<table class='a'>")
        for key in sorted(details['reads'][item]):
            if key in ('problems'):
                continue
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
        HTML.write("</table>\n")
    
    HTML.write("<h3>Pre-Assembly Transformations</h3>\n<div class='a'>\n")
    for line in details["pre-assembly transformation"]:
        HTML.write("<p>"+line+"\n")
    HTML.write("</div>\n")

    if "trim report" in details:
        HTML.write("<h3>Trimming Report</h3>\n<div class='a'>\n")
        for reads in details["trim report"]:
            HTML.write("<b>"+reads+"</b><ul>")
            for report in details["trim report"][reads]:
                HTML.write("<li><a href='%s'>%s</a>\n"%(report, report))
        HTML.write("</div>\n")


    HTML.write("<h3>Assembly</h3>\n")
    HTML.write("<table class='a'>")
    for key in sorted(details['assembly']):
        HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['assembly'][key])))
    HTML.write("</table>\n")

    if "quast_txt" in details:
        HTML.write("<h3>Quast Report</h3>\n")
        HTML.write("<table class='a'>")
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_txt"], "Quast text report"))
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_html"], "Quast html report"))
        HTML.write("</table>\n")
    
    HTML.write("<h3>Post-Assembly Transformations</h3>\n<div class='a'>\n")
    for line in details["post-assembly transformation"]:
        HTML.write("<p>"+line+"<br>\n")
    HTML.write("</div>\n")

    if "Bandage plot" in details:
        HTML.write("<h3>Bandage Plot</h3>\n")
        HTML.write("<div class='a'>")
        HTML.write("<img src='%s'>\n"%(details["Bandage plot"]))
        HTML.write("</div>\n")



def main():
    global Start_time
    Start_time = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='p3_assembly')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', nargs='*', help='list of Pacific Biosciences fastq[.gz] or bam files', required=False, default=[])
    parser.add_argument('--nanopore', nargs='*', help='list of Oxford Nanotech fastq[.gz] or bam files', required=False, default=[])
    parser.add_argument('--sra', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('--recipe', choices=['unicycler', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'auto'], help='assembler to use', default='auto')

    parser.add_argument('--racon_iterations', type=int, default=2, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=2, help='number of times to run pilon per short-read file', required=False)
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--prefix', default='', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', default=Default_genome_size, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', default=300, help='save contigs of this length or longer', required=False)
    parser.add_argument('--min_contig_coverage', default=5, help='save contigs of this coverage or deeper', required=False)
    parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cpus', type=int, default=4)
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
    baseName = "p3_assembly" 
    if len(args.prefix) > 0 and not args.prefix.endswith("_"):
        args.prefix += "_"
    global WorkDir
    WorkDir = baseName+"_work"
    if os.path.exists(WorkDir):
        shutil.rmtree(WorkDir)
    os.mkdir(WorkDir)
    global SaveDir
    SaveDir = os.path.abspath(os.path.join(WorkDir, baseName+"_save"))
    if os.path.exists(SaveDir):
        shutil.rmtree(SaveDir)
    os.mkdir(SaveDir)
    logfileName = os.path.join(SaveDir, args.prefix + baseName + ".log")

    global LOG 
    sys.stderr.write("logging to "+logfileName+"\n")
    LOG = open(logfileName, 'w', 0) #unbuffered 
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(Start_time))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Temporary directory is "+WorkDir+"\n\n")
    LOG.write("Final output will be saved to "+SaveDir+"\n\n")
    details = { 'logfile' : logfileName }
    details["pre-assembly transformation"] = []
    details["post-assembly transformation"] = []
    details["original_items"] = []
    details["reads"] = {}
    details["problem"] = []
    details["platform"] = {}

    if args.illumina:
        platform='illumina'
        for item in args.illumina:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(details, item, platform=platform, interleaved=interleaved)

    if args.iontorrent:
        platform='iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(details, item, platform=platform, interleaved=interleaved)

    if args.pacbio:
        for item in args.pacbio:
            registerReads(details, item, platform='pacbio')

    if args.nanopore:
        for item in args.nanopore:
            registerReads(details, item, platform='nanopore')

    if args.sra:
        fetch_sra_files(args.sra, details)

    if args.anonymous_reads:
        categorize_anonymous_read_files(args.anonymous_reads, details)

    # move into working directory so that all files are local
    original_working_directory = os.getcwd()
    os.chdir(WorkDir)

    if args.trim:
        trimGalore(args, details)

    if args.recipe == "auto":
        #now must decide which assembler to use
        if True:
            # original rule: if any illumina or iontorrent reads present, use Unicycler (long-reads can be present), else use canu for long-reads
            if args.illumina or args.iontorrent:
                args.recipe = "unicycler"
            else:
                args.recipe = "canu"
        else:
            # alternative rule: if any long reads present, use canu
            if args.pacbio or args.nanopore:
                args.recipe = "canu"
            else:
                args.recipe = "unicycler"
    if "spades" in args.recipe:
        contigs = runSpades(args, details)
    elif args.recipe == "unicycler":
        contigs = runUnicycler(args, details)
    elif args.recipe == "canu":
        contigs = runCanu(args, details)
    else:
        LOG.write("cannot interpret args.recipe: "+args.recipe)

    if contigs and os.path.getsize(contigs):
        # now run racon with each long-read file
        for i in range(0, args.racon_iterations):
            for longReadFastq in args.pacbio + args.nanopore:
                LOG.write("runRacon(%s, %s, args, details)\n"%(contigs, longReadFastq))
                raconContigFile = runRacon(contigs, longReadFastq, args, details)
                if raconContigFile is not None:
                    comment = "racon: input %s, output %s"%(contigs, raconContigFile)
                    contigs = raconContigFile
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                else:
                    break
        
    if contigs and os.path.getsize(contigs):
        # now run pilon with each short-read file
        for i in range(0, args.pilon_iterations):
            numChanges = 0
            for shortReadFastq in args.illumina + args.iontorrent:
                LOG.write("runPilon(%s, %s, args, details)\n"%(contigs, shortReadFastq))
                pilonContigFile = runPilon(contigs, shortReadFastq, args, details)
                if pilonContigFile is not None:
                    comment = "pilon: input %s, output %s"%(contigs, pilonContigFile)
                    if os.path.exists(pilonContigFile.replace(".fasta", ".changes")):
                        with open(pilonContigFile.replace(".fasta", ".changes")) as F:
                            lines = F.readlines()
                            numChanges = len(lines)
                            comment += ", num_changes = %d"%numChanges
                    LOG.write(comment+"\n")
                    details["post-assembly transformation"].append(comment)
                    contigs = pilonContigFile
                else:
                    break
            if not numChanges:
                break
        
    if contigs and os.path.getsize(contigs):
        shortReadDepth={}
        longReadDepth={}
        if len(args.illumina + args.iontorrent):
            bamFiles=[]
            for fastq in args.illumina + args.iontorrent:
                bam = runBowtie(contigs, fastq, args.threads, output='bam')
                if bam:
                    bamFiles.append(bam)
            shortReadDepth = calcReadDepth(bamFiles)
        if len(args.pacbio + args.nanopore):
            bamFiles=[]
            for fastq in args.pacbio + args.nanopore:
                bam = runMinimap(contigs, fastq, args.threads, output='bam')
                if bam:
                    bamFiles.append(bam)
            longReadDepth = calcReadDepth(bamFiles)
        saveContigsFile = os.path.join(SaveDir, args.prefix+"contigs.fasta")
        filterContigsByMinLength(contigs, saveContigsFile, args, details, shortReadDepth=shortReadDepth, longReadDepth=longReadDepth)
        runQuast(saveContigsFile, args, details)
        runBandage(args, details)

    htmlFile = os.path.join(SaveDir, args.prefix+"assembly_report.html")
    write_html_report(htmlFile, details)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - Start_time)/3600))
    LOG.close()
    shutil.move(logfileName, os.path.join(args.outputDirectory, os.path.basename(logfileName)))
    fp = file(os.path.join(SaveDir, baseName+".run_details"), "w")
    json.dump(details, fp, indent=2, sort_keys=True)
    fp.close()


if __name__ == "__main__":
    main()
