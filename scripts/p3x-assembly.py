#!/usr/bin/env python
import sys
import subprocess
import argparse
import gzip
import os
import os.path
import re
import shutil
import urllib2
from time import time, localtime, strftime
import json
import sra_tools

"""
This script organizes a command line for either 
Unicycler, or canu as appropriate: 
    canu if only long reads (pacbio or nanopore), 
    Unicycler if illumina or iontorrent 
    and Unicycler for hybrid assemblies, eg Illumina plus PacBio
It can auto-detect read types (illumina, iontorrent, pacbio, nanopore)
It can run Trimmomatic prior to assembling.
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
Path_to_lib = '.'

def parseJasonParameters(args):
    if not os.path.exists(args.params_json):
        raise Exception("cannot find json parameters file %s\n"%args.params_json)
    with open(args.params_json) as json_file:
        data = json.load(json_file)
        args.output_dir = data['output_file']
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

def determineReadFileType(read_id, avgReadLength):
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

def trimPairedReads(readPair, args, illumina=False):
    """ Split files on separator (: or %)
    run Trimmomatic (specifying IlluminaClip if illumina parameter set to true)
    return joined (on separator) paired files, plus one file with unpaired reads """
    LOG.write("trimPairedReads: elapsed seconds = %f\n"%(time()-Start_time))
    m = re.match("^(.+)([%:])(.+)", readPair)
    if not m:
        raise Exception("trimPairedReads cannot split readPair: "+readPair)
    readFile1, separator, readFile2 = m.groups()
    read1_out_base = os.path.basename(readFile1)
    read1_out_base = read1_out_base.split(".")[0]

    read2_out_base = os.path.basename(readFile2)
    read2_out_base = read2_out_base.split(".")[0]

    if read1_out_base == read2_out_base:
        read1_out_base += "_read1"
        read2_out_base += "_read2"
    read1_out_base += "_trim"
    read2_out_base += "_trim"
    read1_out_base = os.path.join(args.output_dir, read1_out_base)
    read2_out_base = os.path.join(args.output_dir, read2_out_base)

    command = ["java", "-jar", os.path.join(Path_to_lib, 'trimmomatic.jar'), "PE", "-threads", str(args.threads)]
    command.extend([readFile1, readFile2, read1_out_base+"_P.fq", read1_out_base+"_U.fq", read2_out_base+"_P.fq", read2_out_base+"_U.fq"])
    if illumina:
        command.append("ILLUMINACLIP:%s:2:30:10"%os.path.join(Path_to_lib, 'illumina_adapters.fa'))
    command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
    command.append("LEADING:%d"%(args.trimmomaticEndQual))
    command.append("TRAILING:%d"%(args.trimmomaticEndQual))
    LOG.write("command = "+" ".join(command)+"\n")
    
    return_code = runProcess(command, shell=False)
    LOG.write("return code = %d\n"%return_code)
    #
    # If we have unpaired reads, copy them to a single file.
    # If we are running in metagenomics mode, discard the unpaired reads since
    # metaSPADES only accepts a single PE library.
    # 
    if os.path.exists(read2_out_base+"_U.fq") and not args.meta:
        UnpairedRead1File = open(read1_out_base+"_U.fq", "a")
        UnpairedRead2File = open(read2_out_base+"_U.fq")
        for line in UnpairedRead2File:
            UnpairedRead1File.write(line)
        UnpairedRead1File.close()
        UnpairedRead2File.close()
    # return value is tuple of two elements
    # first element is the two trimmed paired-read files joined by the separator found on input [: or %]
    # second element is the name of the file with unpaired reads (combining read1 and read2 to one file)
    return(read1_out_base+"_P.fq" + separator + read2_out_base+"_P.fq", read1_out_base+"_U.fq")

def trimSingleReads(readFile, args, illumina=False):
    """ Run Trimmomatic on single read file (include IlluminaClip if illumina==True)"""
    LOG.write("trimSingleReads: elapsed seconds = %f\n"%(time()-Start_time))
    trimmedFileName = os.path.basename(readFile)
    trimmedFileName = trimmedFileName.split('.')[0]
    trimmedFileName += "_trim.fq"
    trimmedFileName = os.path.join(args.output_dir, trimmedFileName)

    command = ["java", "-jar", os.path.join(Path_to_lib, 'trimmomatic.jar'), "SE", "-threads", str(args.threads)]
    command.extend([readFile, trimmedFileName])
    if illumina:
        command.append("ILLUMINACLIP:%s:2:30:10"%os.path.join(Path_to_lib, 'illumina_adapters.fa'))
    command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
    command.append("LEADING:%d"%(args.trimmomaticEndQual))
    command.append("TRAILING:%d"%(args.trimmomaticEndQual))
    LOG.write("command = "+" ".join(command)+"\n")
    return_code = runProcess(command, shell=False)
    LOG.write("return code = %d\n"%return_code)
    return trimmedFileName

def verifyReadPairing(readPair, output_dir):
    """ Read both files, write 3 new files: paired1, paired2, unpaired
    return paired1[:%]paired2, unpaired
    """
    m = re.match("^(.+)([%:])(.+)", readPair)
    if not m:
        raise Exception("verifyReadPairing cannot split readPair: "+readPair)
    readFile1, separator, readFile2 = m.groups()
    if readFile1.lower().endswith(".gz"):
        F = gzip.open(readFile1)
    else:
        F = open(readFile1)
    i = 0
    id = None
    read1 = {}
    seqqual = '' # for 3-line record of sequence and quality
    for line in F:
        if i % 4 == 0:
            if id:
                read1[id] = seqqual
            seqqual = ''
            id = line.split()[0]
        else:
            seqqual += line
        i += 1
    F.close()
    if id:
        read1[id] = seqqual

    basename = os.path.join(output_dir, os.path.basename(readFile1))
    verifiedPairedFile1 = basename + "_paired.fq"
    unpairedReadFile = basename + "_reads12_unpaired.fq"
    basename = os.path.join(output_dir, os.path.basename(readFile2))
    verifiedPairedFile2 = basename + "_paired.fq"
    PairedOut1 = open(verifiedPairedFile1, 'w')
    PairedOut2 = open(verifiedPairedFile2, 'w')
    UnpairedOut = open(unpairedReadFile, 'w')

    if readFile2.lower().endswith(".gz"):
        F = gzip.open(readFile2)
    else:
        F = open(readFile2)
    i = 0
    found = False
    for line in F:
        if i % 4 == 0:
            id = line.split()[0]
            found = id in read1
            if not found:
                idmod = re.sub("2$", "1", id)
                found = idmod in read1
                if found:
                    id = idmod
            if not found:
                idmod = re.sub("2:(\d+)$", "1:\1", id)
                found = idmod in read1
                if found:
                    id = idmod
            if not found:
                idmod = re.sub("1$", "2", id)
                found = idmod in read1
                if found:
                    id = idmod
            if not found:
                idmod = re.sub("1:(\d+)$", "2:\1", id)
                found = idmod in read1
                if found:
                    id = idmod
            if found:
                PairedOut1.write(id+"\n"+read1[id])
                PairedOut2.write(id+"\n")
                del read1[id]
            else:
                UnpairedOut.write(id+"\n")
        else:
            if found:
                PairedOut2.write(line)
            else:
                UnpairedOut.write(line)
        i += 1
    PairedOut1.close()
    PairedOut2.close()

    for id in read1:
        UnpairedOut.write(id+"\n"+read1[id])
    UnpairedOut.close()
    # separator is ':' or '%' from input
    return (verifiedPairedFile1+separator+verifiedPairedFile2, unpairedReadFile)

def studyReadFile(filename, details=None):
    LOG.write("studyReadFile(%s): elapsed seconds = %f\n"%(filename, time()-Start_time))
    #return(Read_file_type[filename], Read_id_sample[filename], Avg_read_length[filename])
    # figures out Read_file_type, collects a sample of read IDs, and average read length
    Read_file_type[filename] = 'na'
    Read_id_sample[filename] = []
    if filename.endswith("gz"):
        F = gzip.open(filename)
    else:
        F = open(filename)
    text = F.read(Default_bytes_to_sample) #read X number of bytes for text sample
    F.close()
    LOG.write("  file text sample %s:\n%s\n\n"%(filename, text[0:50]))
    Read_id_sample[filename] = []
    lines = text.split("\n")
    if len(lines) < 2:
        raise Exception("text sample (length %d) lacks at least 2 lines"%len(text))
    read_format = None
    if lines[0].startswith("@"):
        read_format = 'fastq'
    elif lines[0].startswith(">"):
        read_format = 'fasta'
    Avg_read_length[filename] = 0
    readLengths = []
    if read_format == 'fastq':
        for i, line in enumerate(lines):
            if i % 4 == 0:
                Read_id_sample[filename].append(line.split(' ')[0]) # get part up to first space, if any 
            elif i % 4 == 1:
                readLengths.append(len(line)-1)
    elif read_format == 'fasta':
        read = ''
        for line in lines:
            if line.startswith(">"):
                Read_id_sample[filename].append(line[1:].split(' ')[0]) # get part after '>' and up to first space, if any 
                if len(read) > 0:
                    readLengths.append(len(read))
                read = ''
            else:
                read += line.rstrip()
    Avg_read_length[filename] = sum(readLengths)/float(len(readLengths))
    Read_file_type[filename] = determineReadFileType(lines[0], Avg_read_length[filename])
    LOG.write("found read type %s average read length %.1f\n"%(Read_file_type[filename], Avg_read_length[filename]))
    if details:
        if not "inferred_read_type" in details:
            details["inferred_read_type"] = {}
        details["inferred_read_type"][filename] = Read_file_type[filename]
        if not "estimated_read_length" in details:
            details["estimated_read_length"] = {}
        details["estimated_read_length"][filename] = Avg_read_length[filename]
    return(Read_file_type[filename], Read_id_sample[filename], Avg_read_length[filename])

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
    LOG.write("categorize_anonymous_read_files: elapsed seconds = %f\n"%(time()-Start_time))
    LOG.write("  files=%s\n"%("\t".join(args.anonymous_reads)))
    for filename in args.anonymous_reads:
        studyReadFile(filename, details)

    # try to find paired files
    membersOfPairs = set()
    for i, filename1 in enumerate(args.anonymous_reads[:-1]):
        for filename2 in args.anonymous_reads[i+1:]:
            charDiffs = findSingleDifference(filename1, filename2)
# charDiffs will be not None if the strings match at all but one character (presumably '1' vs '2')
            if charDiffs:
                try:
                    intDiffs = (int(float(charDiffs[0])), int(float(charDiffs[1])))
                    pairedFiles = None
                    if intDiffs[0] == 1 and intDiffs[1] == 2:
                        pairedFiles = (filename1, filename2)
                    elif intDiffs[1] == 1 and intDiffs[0] == 2:
                        pairedFiles = (filename2, filename1)
                    if pairedFiles:
                        if Read_file_type[filename1] != Read_file_type[filename2]:
                            LOG.write("!Discordant fileTypes for %s vs %s, %s vs %s\n"%(filename1, filename2, Read_file_type[filename1], Read_file_type[filename2]))
                        if Read_file_type[filename1] == 'illumina':
                            if not args.illumina:
                                args.illumina = []
                            args.illumina.append(":".join(pairedFiles))
                            LOG.write("appending to args.illumina: %s %s\n"%pairedFiles)
                        elif Read_file_type[filename1] == 'iontorrent':
                            if not args.iontorrent:
                                args.iontorrent = []
                            args.iontorrent.append(":".join(pairedFiles))
                            LOG.write("appending to args.iontorrent: %s %s\n"%pairedFiles)
                        else: # neither illumina vs iontorrent, perhaps 'sra'
                            if Avg_read_length[filename1] < Max_short_read_length:
                                #call it illumina
                                if not args.illumina:
                                    args.illumina=[]
                                args.illumina.append(":".join(pairedFiles))
                                LOG.write("Calling file pair %s %s, mean length %d, to be 'illumina', from %s\n"%(filename1, filename2, Avg_read_length[filename1], Read_file_type[filename1]))
                        membersOfPairs.add(pairedFiles[0])
                        membersOfPairs.add(pairedFiles[1])
                except:
                    pass

    for filename in args.anonymous_reads:
        if filename not in membersOfPairs:
            if Read_file_type[filename] == 'illumina':
                if not args.illumina:
                    args.illumina = []
                args.illumina.append(filename)
                LOG.write("appending to args.illumina: %s\n"%filename)
            elif Read_file_type[filename] == 'iontorrent':
                if not args.iontorrent:
                    args.iontorrent = []
                args.iontorrent.append(filename)
                LOG.write("appending to args.iontorrent: %s\n"%filename)
            elif Read_file_type[filename] == 'pacbio':
                if not args.pacbio:
                    args.pacbio = []
                args.pacbio.append(filename)
                LOG.write("appending to args.pacbio: %s\n"%filename)
            elif Read_file_type[filename] == 'nanopore':
                if not args.nanopore:
                    args.nanopore = []
                args.nanopore.append(filename)
                LOG.write("appending to args.nanopore: %s\n"%filename)
            elif Avg_read_length[filename] < Max_short_read_length:
                #call it illumina
                if not args.illumina:
                    args.illumina=[]
                args.illumina.append(filename)
                LOG.write("Calling file %s, mean length %d, to be 'illumina', from %s\n"%(filename, Avg_read_length[filename], Read_file_type[filename]))
            else:
                if not args.pacbio:
                    args.pacbio=[]
                args.pacbio.append(filename)
                LOG.write("Calling file %s, mean length %d, to be 'pacbio', from %s\n"%(filename, Avg_read_length[filename], Read_file_type[filename]))
    return

def fetch_sra_files(args, details):
    """ Need to change this to call external library Zane and Andew wrote
    Use ftp to get all SRA files.
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    Append to appropriate parts of args (e.g., args.illumina or args.iontorrent).
    """
    LOG.write("fetch_sra_files: elapsed seconds = %f\n"%(time()-Start_time))
    LOG.write("fetch_sra_files() args.sra="+" ".join(args.sra)+"\n")
    if not args.sra:
        return
    for sra in args.sra:
        sra_dir = ""
        if sra.endswith(".sra"):
            sra_dir, sra = os.path.split(sra)
            sra = sra[:-4] # trim off trailing ".sra", will later detect presence of "xxx.sra" file if it exists

        runinfo = sra_tools.get_runinfo(sra)
        if runinfo['Platform'] == 'PACBIO_SMRT':
            if not args.pacbio:
                args.pacbio = []
            listToAddTo = args.pacbio
        elif runinfo['Platform'] == "ILLUMINA":
            if not args.illumina:
                args.illumina = []
            listToAddTo = args.illumina
        elif runinfo['Platform'] == "ION_TORRENT":
            if not args.iontorrent:
                args.iontorrent = []
            listToAddTo = args.iontorrent
        elif runinfo['Platform'] == "OXFORD_NANOPORE":
            if not args.nanopore:
                args.nanopore = []
            listToAddTo = args.nanopore
        else:
            raise Exception("Problem: Cannot process sra data from platform %s\n"%runinfo['Platform'])

        if not os.path.exists(os.path.join(sra_dir, sra+".sra")):
            LOG.write("downloading %s\n"%sra)
            sra_tools.ftp_download_single_run(sra)

        if not os.path.exists(os.path.join(sra_dir, sra+".sra")):
            raise Exception("Problem: file %s.sra does not exist after trying to download %s\n"%(sra, sra_file_url))

        if runinfo['LibraryLayout'].startswith("SINGLE"):
            sra_tools.fastqDumpExistingSraFile(os.path.join(sra_dir, sra+".sra"), splitFiles=False)
            #runProcess(["fastq-dump", sra+".sra"], shell=False)
            if not os.path.exists(sra+".fastq"):
                raise Exception("Problem: file %s.fastq does not exist after running fastq-dump on %s.sra\n"%(sra, sra))
            listToAddTo.append(sra+".fastq")
        elif runinfo['LibraryLayout'].startswith("PAIRED"):
            sra_tools.fastqDumpExistingSraFile(os.path.join(sra_dir, sra+".sra"), splitFiles=True)
            #runProcess(["fastq-dump", "--split-files", sra+".sra"], shell=False)
            if not (os.path.exists(sra+"_1.fastq") and os.path.exists(sra+"_2.fastq")):
                raise Exception("Problem: file %s_1.fastq and/or %s_2.fastq do not exist after running fastq-dump --split-files on %s.sra\n"%(sra, sra, sra))
            listToAddTo.append(sra + "_1.fastq:" + sra +"_2.fastq")
    return

def study_all_read_files(args, details):
    LOG.write("study_all_read_files: elapsed seconds = %f\n"%(time()-Start_time))
    fileItemType={}
    filePairs=[]
    singleFiles=[]
    if args.illumina:
        for item in args.illumina:
            fileItemType[item] = 'illumina'
            if ':' in item or '%' in item:
                filePairs.append(item)
            else:
                singleFiles.append(item)
    if args.iontorrent:
        for item in args.iontorrent:
            fileItemType[item] = 'iontorrent'
            if ':' in item or '%' in item:
                filePairs.append(item)
            else:
                singleFiles.append(item)
    if args.pacbio:
        for item in args.pacbio:
            fileItemType[item] = 'pacbio'
            singleFiles.append(item)
    if args.nanopore:
        for item in args.nanopore:
            fileItemType[item] = 'nanopore'
            singleFiles.append(item)

    for item in filePairs:
        LOG.write("studying read pair %s\n"%item)
        LOG.write("\tclaimed type = %s\n"%fileItemType[item])
        pair = item.split(':')
        if '%' in item:
            pair = item.split('%')
        charDiffs = findSingleDifference(pair[0], pair[1])
        if charDiffs:
            LOG.write("\tsingle char diff=%s %s\n"%(charDiffs[0], charDiffs[1]))
        else:
            LOG.write("\tno single char difference found\n")
        if pair[0] not in Read_file_type:
            studyReadFile(pair[0], details)
        if pair[1] not in Read_file_type:
            studyReadFile(pair[1], details)
        if Read_file_type[pair[0]] == Read_file_type[pair[1]]:
            LOG.write("\tinspected file types congruent: %s\n"%Read_file_type[pair[0]])
            if Read_file_type[pair[0]] != fileItemType[item]:
                LOG.write("\t!discrepancy with claimed type: %s is type %s\n"%(pair[1], Read_file_type[pair[1]]))
        else:
            LOG.write("\t!inspected file types incongurent: %s vs %s\n"%(Read_file_type[pair[0]], Read_file_type[pair[1]]))
        allPairsMatch = True
        for read_idPair in zip(Read_id_sample[pair[0]], Read_id_sample[pair[1]]):
            if read_idPair[0] != read_idPair[1]:
                allPairsMatch = False
                
        LOG.write("\tread IDs tested for match for files %s "%str(pair)+" result = %s\n"%str(allPairsMatch))

    for filename in singleFiles:
        studyReadFile(filename, details)
        if Read_file_type[filename] != fileItemType[filename]:
            LOG.write("\t!discrepancy with claimed type: %s is type %s\n"%(filename, Read_file_type[filename]))
    return

def writeSpadesYamlFile(args):
    LOG.write("writeSpadesYamlFile: elapsed seconds = %f\n"%(time()-Start_time))
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    outfileName = os.path.join(args.output_dir, "spades_yaml_file.txt")
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    LOG.write("illumina: "+", ".join(args.illumina)+"\n")
    LOG.write("iontorrent: "+", ".join(args.iontorrent)+"\n")
    
    single_end_reads = []
    paired_end_reads = [[], []]
    mate_pair_reads = [[], []]
    all_read_files = []
    if args.illumina:
        all_read_files = args.illumina
    elif args.iontorrent:
        all_read_files = args.iontorrent
    print all_read_files
    for item in all_read_files:
        if ":" in item:
            pair = item.split(":")
            if os.path.exists(pair[0]) and os.path.getsize(pair[0]) > 0:
                f = os.path.abspath(pair[0])
                paired_end_reads[0].append(f)
                f = os.path.abspath(pair[1])
                paired_end_reads[1].append(f)
        elif "%" in item:
            pair = item.split("%")
            if os.path.exists(pair[0]) and os.path.getsize(pair[0]) > 0:
                f = os.path.abspath(pair[0])
                mate_pair_reads[0].append(f)
                f = os.path.abspath(pair[1])
                mate_pair_reads[1].append(f)
        else:
            if os.path.exists(item) and os.path.getsize(item) > 0:
                single_end_reads.append(os.path.abspath(item))


    #
    # Validate arguments for metagenomic spades. It can only run with
    # a single paired-end library.
    #

    if args.meta:
        if len(single_end_reads) > 0 or len(paired_end_reads[0]) > 1:
            sys.stderr.write("SPAdes in metagenomics mode can only process a single paired-end read file\n")
            sys.exit(1);

    precedingElement=False
    if len(single_end_reads):
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(single_end_reads))
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
    if args.pacbio:
        if precedingElement:
            OUT.write(",\n")
        pacbio_reads = []
        for f in args.pacbio:
            pacbio_reads.append(os.path.abspath(f))
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True
    if args.nanopore:
        if precedingElement:
            OUT.write(",\n")
        nanopore_reads = []
        for f in args.nanopore:
            nanopore_reads.append(os.path.abspath(f))
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement=True

    OUT.write("]\n")
    OUT.close()
    return(outfileName)    

def filterContigsByMinLength(inputFile, outputFile, min_length = 0):
    """ Write only sequences at or above min_length to output file."""
    LOG.write("filterContigsByMinLength(%s, %s, %d) elapsed seconds = %f\n"%(inputFile, outputFile, min_length, time()-Start_time))
    with open(inputFile) as IN:
        with open(outputFile, 'w') as OUT:
            seqId=None
            seq = ""
            for line in IN:
                if line.startswith(">"):
                    if seqId and len(seq) >= min_length:
                        OUT.write(">%s\n%s\n"%(seqId, seq))
                        seqId = None
                        seq = ""
                    m = re.match(">(\S+)", line)
                    seqId = m.group(1)
                else:
                    seq += line.rstrip()
            if seqId and len(seq) >= min_length:
                OUT.write(">%s\n%s\n"%(seqId, seq))

def runUnicycler(args, details):
    LOG.write("runUnicycler: elapsed seconds = %f\n"%(time()-Start_time))
    command = ["unicycler", "-t", str(args.threads), "-o", args.output_dir, "--pilon_path", args.pilon_path]
    if args.min_contig_length:
        command.extend(("--min_fasta_length", str(args.min_contig_length)))
    command.extend(("--keep", "0")) # keep only assembly.gfa, assembly.fasta and unicycler.log
    if args.illumina:
        for item in args.illumina:
            if ":" in item:
                read1, read2 = item.split(":")
                command.extend(("-1", read1, "-2", read2))
            else:
                command.extend(("-s", item))
    # it is not quite right to send iontorrent data to spades through unicycler because the --iontorrent flag to spades will not be set
    if args.iontorrent:
        for item in args.iontorrent:
            if ":" in item:
                read1, read2 = item.split(":")
                command.extend(("-1", read1, "-2", read2))
            else:
                command.extend(("-s", item))
    if args.pacbio:
        for item in args.pacbio:
            command.extend(("-l", item))
    if args.nanopore:
        for item in args.nanopore:
            command.extend(("-l", item))

    LOG.write("Unicycler command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    unicyclerStartTime = time()

    return_code = runProcess(command, shell=False)
    LOG.write("return code = %d\n"%return_code)

    unicyclerEndTime = time()
    elapsedTime = unicyclerEndTime - unicyclerStartTime

    details.update( { 'start_time': unicyclerStartTime,
                'end_time': unicyclerEndTime,
                'elapsed_time' : elapsedTime,
                'assembler': 'unicycler',
                'command_line': command,
                'output_path': os.path.abspath(args.output_dir)
                })
    LOG.write("Duration of Unicycler run was %.1f hours\n"%(elapsedTime/3600.0))
    unicyclerLogFile = args.prefix+"unicycler.log"
    shutil.move(os.path.join(args.output_dir, "unicycler.log"), os.path.join(args.output_dir, unicyclerLogFile))
    details['output_files'] = []
    details['output_files'].append([unicyclerLogFile, 'txt', 'Unicycler logfile'])

    assemblyFile = args.prefix+"assembly.fasta"
    assemblyGraphFile = args.prefix+"assembly_graph.gfa"
    filterContigsByMinLength(os.path.join(args.output_dir, "assembly.fasta"), os.path.join(args.output_dir, assemblyFile), args.min_contig_length)
    shutil.move(os.path.join(args.output_dir, "assembly.gfa"), os.path.join(args.output_dir, assemblyGraphFile))
    details['output_files'].append([assemblyFile, 'fasta', 'Genome assembly'])
    details['output_files'].append([assemblyGraphFile, 'gfa', 'Assembly graph'])
    
    if not args.no_quast:
        quastDir = os.path.join(args.output_dir, "quast_out")
        quastCommand = [args.quast_exec,
                        "-o", quastDir,
                        "-t", str(args.threads),
                        "--min-contig", str(args.min_contig_length),
                        "--gene-finding",
                        os.path.join(args.output_dir, assemblyFile)]
        LOG.write("running quast: "+" ".join(quastCommand)+"\n")
        return_code = runProcess(quastCommand, shell=False)
        LOG.write("return code = %d\n"%return_code)
        shutil.move(os.path.join(quastDir, "report.html"), os.path.join(args.output_dir, args.prefix+"quast_report.html"))
        shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(args.output_dir, args.prefix+"quast_report.txt"))
        shutil.move(os.path.join(quastDir, "icarus_viewers"), os.path.join(args.output_dir, "icarus_viewers"))
        details['output_files'].append([args.prefix+"quast_report.html", 'html', 'Quast report'])
        details['output_files'].append([args.prefix+"quast_report.txt", 'txt', 'Quast report'])
        details['output_files'].append([args.prefix+"icarus_viewers", 'dir', 'Quast icarus viewers'])

    LOG.write("Duration of Unicycler run was %f hours\n"%(elapsedTime/3600.0))

def runSpades(args, details):
    LOG.write("runSpades: elapsed seconds = %f\n"%(time()-Start_time))
    #if ("illumina_pe" in args.output_dirr "illumina_se" in args) and ("iontorrent_pe" in args or "iontorrent_se" in args):
    if args.illumina and args.iontorrent:
        details["Fatal_error"] = "SPAdes cannot process both Illumina and IonTorrent reads in the same run"
        return
    command = ["spades.py", "--threads", str(args.threads), "-o", args.output_dir]
    if args.singlecell:
        command.append("--sc")
    if args.iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(args)
    command.extend(["--dataset", yamlFile])
    if args.only_assembler:
        command.append("--only-assembler")
    if args.trusted_contigs:
        command.extend(["--trusted-contigs", args.trusted_contigs])
    if args.untrusted_contigs:
        command.extend(["--untrusted-contigs", args.untrusted_contigs])
    if not args.no_careful and not args.meta:
        command.append("--careful")
    if args.memory:
        command.extend(["-m", str(args.memory)])
    if args.meta:
        command.append("--meta")
    LOG.write("SPAdes command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    spadesStartTime = time()

    return_code = runProcess(command, shell=False)
    LOG.write("return code = %d\n"%return_code)

    spadesEndTime = time()
    elapsedTime = spadesEndTime - spadesStartTime

    details.update( { 'spades_start_time': spadesStartTime,
                'spades_end_time': spadesEndTime,
                'spades_elapsed_time' : elapsedTime,
                'assembler': 'spades',
                'command_line': command,
                'output_path': os.path.abspath(args.output_dir)
                })
    scaffoldsFile = args.prefix+"scaffolds.fasta"
    contigsFile = args.prefix+"contigs.fasta"
    assemblyGraphFile = args.prefix+"assembly_graph.gfa"
    spadesLogFile = args.prefix+"spades.log"
    
    shutil.move(os.path.join(args.output_dir, "spades.log"), os.path.join(args.output_dir, spadesLogFile))
    details['output_files'] = []
    details['output_files'].append([spadesLogFile, 'txt', 'Spades logfile'])
    
    if not os.path.exists(os.path.join(args.output_dir, "contigs.fasta")):
        LOG.write("Duration of SPAdes run was %.1f hours\n"%(elapsedTime/3600.0))
        LOG.write("spades failed to generate contigs.fasta")
        return
    filterContigsByMinLength(os.path.join(args.output_dir, "contigs.fasta"), os.path.join(args.output_dir, contigsFile), args.min_contig_length)
    filterContigsByMinLength(os.path.join(args.output_dir, "scaffolds.fasta"), os.path.join(args.output_dir, scaffoldsFile), args.min_contig_length)
    shutil.move(os.path.join(args.output_dir, "assembly_graph_with_scaffolds.gfa"), os.path.join(args.output_dir, assemblyGraphFile))
    details['output_files'].append([scaffoldsFile, 'fasta', 'Generated scaffolds'])
    details['output_files'].append([contigsFile, 'fasta', 'Generated contigs'])
    details['output_files'].append([assemblyGraphFile, 'gfa', 'Assembly graph'])
    
    if not args.no_quast:
        quastDir = os.path.join(args.output_dir, "quast_out")
        quastCommand = [args.quast_exec,
                        "-o", quastDir,
                        "-t", str(args.threads),
                        "--min-contig", str(args.min_contig_length),
                        "--gene-finding",
                        os.path.join(args.output_dir, scaffoldsFile),
                        os.path.join(args.output_dir, contigsFile)]
        LOG.write("running quast: "+" ".join(quastCommand)+"\n")
        return_code = runProcess(quastCommand, shell=False)
        LOG.write("return code = %d\n"%return_code)
        shutil.move(os.path.join(quastDir, "report.html"), os.path.join(args.output_dir, args.prefix+"quast_report.html"))
        shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(args.output_dir, args.prefix+"quast_report.txt"))
        shutil.move(os.path.join(quastDir, "icarus_viewers"), os.path.join(args.output_dir, "icarus_viewers"))
        details['output_files'].append([args.prefix+"quast_report.html", 'html', 'Quast report'])
        details['output_files'].append([args.prefix+"quast_report.txt", 'txt', 'Quast report'])
        details['output_files'].append([args.prefix+"icarus_viewers", 'dir', 'Quast icarus viewers'])

    LOG.write("Duration of SPAdes run was %.1f hours\n"%(elapsedTime/3600.0))
    return

def runCanu(args, details):
    LOG.write("runCanu: elapsed seconds = %f\n"%(time()-Start_time))
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
    command = ["canu", "-d", args.output_dir, "-p", "canu", "useGrid=false", "genomeSize=%s"%args.genome_size]
    command.extend(["maxMemory=" + str(args.memory), "maxThreads=" + str(args.threads)])
    """
    The 'minMemory', 'maxMemory', 'minThreads' and 'maxThreads'
    options will apply to all jobs, and can be used to artificially
    limit canu to a portion of the current machine. 
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """

#    for prefix in ("mhap", "mmap", "ovl", "ovb", "cor", "red", "oea", "bat",
#            "cns"):
#        command.append(prefix+"Concurrency=1")
    if args.pacbio:
        command.append("-pacbio-raw")
        command.extend(args.pacbio) #allow multiple files
    if args.nanopore:
        command.append("-nanopore-raw")
        command.extend(args.nanopore) #allow multiple files
    LOG.write("canu command =\n"+" ".join(command)+"\n")
    #LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")

    canuStartTime = time()
    canuLogFile = open(os.path.join(args.output_dir, "canu.log"), "w")
    return_code = runProcess(command, shell=False, stderr=canuLogFile)
    LOG.write("return code = %d\n"%return_code)
    canuEndTime = time()
    elapsedTime = canuEndTime - canuStartTime

    details.update( { 'canu_start_time': canuStartTime,
                'canu_end_time': canuEndTime,
                'canu_elapsed_time' : elapsedTime,
                'assembler': 'canu',
                'command_line': command,
                'output_path': os.path.abspath(args.output_dir)
                })
    
    details['output_files'] = []
    if os.path.exists(os.path.join(args.output_dir, "canu.log")):
        details['output_files'].append(["canu.log", 'txt', 'Canu log'])
        if os.path.exists(os.path.join(args.output_dir, "canu.report")):
            canuReportFile = args.prefix+"canu_report.txt"
            shutil.move(os.path.join(args.output_dir, "canu.report"), os.path.join(args.output_dir, canuReportFile))
            details['output_files'].append([canuReportFile, 'txt', 'Canu report'])

    if not os.path.exists(os.path.join(args.output_dir, "canu.contigs.fasta")):
        LOG.write("Duration of canu run was %.1f hours\n"%(elapsedTime/3600.0))
        LOG.write("canu failed to generate contigs.fasta\n")
        return
    unitigsFile = args.prefix+"unitigs.fasta"
    contigsFile = args.prefix+"contigs.fasta"
    unitigsGraphFile = args.prefix+"unitigs.gfa"
    contigsGraphFile = args.prefix+"contigs.gfa"
    filterContigsByMinLength(os.path.join(args.output_dir, "canu.contigs.fasta"), os.path.join(args.output_dir, contigsFile), args.min_contig_length)
    filterContigsByMinLength(os.path.join(args.output_dir, "canu.unitigs.fasta"), os.path.join(args.output_dir, unitigsFile), args.min_contig_length)
    shutil.move(os.path.join(args.output_dir, "canu.unitigs.gfa"), os.path.join(args.output_dir, unitigsGraphFile))
    shutil.move(os.path.join(args.output_dir, "canu.contigs.gfa"), os.path.join(args.output_dir, contigsGraphFile))
    details['output_files'].append([unitigsFile, 'fasta', 'Canu unitigs'])
    details['output_files'].append([contigsFile, 'fasta', 'Canu contigs'])
    details['output_files'].append([unitigsGraphFile, 'gfa', 'Unitigs graph'])
    details['output_files'].append([contigsGraphFile, 'gfa', 'Contigs graph'])
    if not args.no_quast:
        quastDir = os.path.join(args.output_dir, "quast_out")
        quastCommand = [args.quast_exec,
                        "-o", quastDir,
                        "-t", str(args.threads),
                        "--min-contig", str(args.min_contig_length),
                        "--gene-finding",
                        os.path.join(args.output_dir, unitigsFile),
                        os.path.join(args.output_dir, contigsFile)]
        LOG.write("running quast: "+" ".join(quastCommand)+"\n")
        return_code = runProcess(quastCommand, shell=False)
        LOG.write("return code = %d\n"%return_code)
        shutil.move(os.path.join(quastDir, "report.html"), os.path.join(args.output_dir, args.prefix+"quast_report.html"))
        shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(args.output_dir, args.prefix+"quast_report.txt"))
        shutil.move(os.path.join(quastDir, "icarus_viewers"), os.path.join(args.output_dir, "icarus_viewers"))
        details['output_files'].append([args.prefix+"quast_report.html", 'html', 'Quast report'])
        details['output_files'].append([args.prefix+"quast_report.txt", 'txt', 'Quast report'])
        details['output_files'].append([args.prefix+"icarus_viewers", 'dir', 'Quast icarus viewers'])

    LOG.write("Duration of canu run was %.1f hours\n"%(elapsedTime/3600.0))

def main():
    global Path_to_lib
    #Path_to_lib = os.path.dirname(sys.argv[0])
    #Path_to_lib = "/".join(Path_to_lib.split("/")[:-1])
    #Path_to_lib += "/lib"
    global Start_time
    Start_time = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output_dir', default='.', help='output directory.', required=False)
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or "%%" between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', nargs='*', help='list of IonTorrent[.gz] files or pairs, ":" between paired-end-files', required=False, default=[])
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--meta', action = 'store_true', help='run in metagenome mode', required=False)
    parser.add_argument('--pacbio', nargs='*', help='list of Pacific Biosciences fastq[.gz] or bam files', required=False)
    parser.add_argument('--sra', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--nanopore', nargs='*', help='list of Oxford Nanotech fastq[.gz] or bam files', required=False)
    parser.add_argument('--prefix', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', default=Default_genome_size, help='expected genome size for canu, e.g. 300k, 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', default=200, type=int, help='save contigs of this length or longer', required=False)
    #parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--anonymous_reads', nargs='*', help="unspecified read files, types automatically inferred.")
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('--no_careful', action = 'store_true', help='turn off careful flag to SPAdes (faster)', required=False)
    parser.add_argument('--spades_only', action = 'store_true', help='Use SPAdes instead of Unicycler for short reads or hybrid', required=False)
    parser.add_argument("--only-assembler", action = 'store_true', help='turn off SPAdes error correction step', required=False)
    parser.add_argument('-t', '--threads', metavar='cpus', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, help='RAM limit for SPAdes in Gb', default=250)
    parser.add_argument('--bytes_to_sample', metavar='bytes', type=int, default=Default_bytes_to_sample, help='how much to sample from read files to test file type')
    parser.add_argument('--runTrimmomatic', action = 'store_true', help='run trimmomatic on Illumina or Iontorrent fastq files')
    #parser.add_argument('--trimmomatic_jar', default='trimmomatic.jar', help='trimmomatic jar file, with path')
    parser.add_argument('--illuminaAdapters', default='illumina_adapters.fa', help='illumina adapters file, looked for in script directory')
    parser.add_argument('--trimmomaticWindow', metavar="width", type=int, default=Default_window_width, help='window width for trimming')
    parser.add_argument('--trimmomaticMinQual', metavar="phred", type=int, default=Default_window_quality, help='score of window below which 3\' end is trimmed')
    parser.add_argument('--trimmomaticEndQual', metavar="phred", type=int, default=Default_end_quality, help='score at which individual 3\' bases are trimmed')
    parser.add_argument('--no_quast', action = 'store_true', help='turn off runing quast for assembly quality statistics')
    parser.add_argument('--quast_exec', default='quast.py', help='path to quast.py (if not on path)')
    parser.add_argument('--pilon_path', default='pilon', help='path to pilon executable or jar')
    parser.add_argument('--run-details', help='JSON-format document describing details of the run', required=False)
    parser.add_argument('--logfile', help='Log file', required=False)
    parser.add_argument('--params_json', help="JSON file with additional information.")
    parser.add_argument('--path', help="Add the given directories to the PATH", nargs='*', required=False)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    if args.params_json:
        parseJsonParameters(args)
    args.prefix = os.path.basename(args.output_dir)
    args.prefix = args.prefix.rstrip("_")+"_"
    if args.logfile:
        logfileName = os.path.abspath(args.logfile)
    else:
        logfileName = os.path.basename(sys.argv[0])
        logfileName = logfileName.replace(".py", "")
        logfileName = args.prefix + logfileName
        logfileName = os.path.join(args.output_dir, logfileName)+".log"
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    global LOG 
    sys.stderr.write("logging to "+logfileName+"\n")
    LOG = open(logfileName, 'w', 0) #unbuffered 
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(Start_time))+"\n")
    LOG.write("args= "+str(args)+"\n\n")

    if args.path:
        os.environ["PATH"] = ":".join(args.path) + ":" + os.environ["PATH"]

    details = { 'logfile' : logfileName }
    if args.anonymous_reads:
        categorize_anonymous_read_files(args, details)
# if any illumina or iontorrent reads present, use Unicycler (long-reads can be present), else use canu for long-reads
    if args.sra:
        fetch_sra_files(args, details)
    study_all_read_files(args, details)
    if args.illumina or args.iontorrent:
        for fileList in (args.illumina, args.iontorrent):
            if not fileList:
                continue
            processedFileList = []
            for item in fileList:
                if ':' in item or '%' in item:
                    (verifiedPair, unpaired) = verifyReadPairing(item, args.output_dir)
                    if args.runTrimmomatic:
                        (trimmedPair, trimmedUnpaired) = trimPairedReads(verifiedPair, args, illumina=(fileList == args.illumina))
                        processedFileList.append(trimmedPair)
                        #
                        # For metagenomics spades, we can't use the unpaired reads.
                        #
                        if not args.meta:
                            processedFileList.append(trimmedUnpaired)
                        trimmedUnpaired = trimSingleReads(unpaired, args, illumina=(fileList == args.illumina))
                    else:
                        processedFileList.append(verifiedPair)
                        if os.stat(unpaired).st_size > 0:
                            processedFileList.append(unpaired)
                else:
                    if args.runTrimmomatic:
                        trimmedSingleReads = trimSingleReads(item, args, illumina=(fileList == args.illumina))
                        processedFileList.append(trimmedSingleReads)
                    else:
                        processedFileList.append(item)
            if fileList == args.illumina:
                args.illumina = processedFileList
            else:
                args.iontorrent = processedFileList
        if args.spades_only:
            runSpades(args, details)
        else:
            runUnicycler(args, details)
    else: # long reads only, use canu
        runCanu(args, details)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    if args.run_details:
        LOG.write("writing details in json format to %s\n"%args.run_details)
        if "output_files" not in details:
            details["output_files"] = []
        details["output_files"].append([args.run_details, "json", "run_details"])
        fp = file(os.path.join(args.output_dir, args.run_details), "w")
        json.dump(details, fp, indent=2)
        fp.close()

def fatalError(msg):
    LOG.write("Fatal error: %s\n" % msg)
    sys.stderr.write("Fatal error: %s\n" % msg)
    sys.exit(1)

def runProcess(cmd, **opts):
    try:
        rc = subprocess.call(cmd, **opts)
        return rc
    except OSError as e:
        fatalError("Error running %s: %s" % (cmd, e))
    except Exception as e:
        fatalError("Error running %s: %s" % (cmd, e))

if __name__ == "__main__":
    main()
 
