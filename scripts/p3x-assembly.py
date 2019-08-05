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

def determineReadType(read_id, avgReadLength):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    LOG.write("determineReadType() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(Start_time))+"\n")
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
    LOG.write("trimGalore() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    for platformReads in (args.illumina, args.iontorrent):
        for index, reads in enumerate(platformReads): 
            command = ['trim_galore', '-j', str(args.threads), '-o', WorkDir]
            if ':' in reads:
                splitReads = reads.split(":")
                #-j 4 -o testTrim --length 30 --paired SRR1395326_1_10pct.fastq SRR1395326_2_10pct.fastq
                command.extend(["--paired", reads[0], reads[1]])
                LOG.write("command: "+" ".join(command))
                return_code = subprocess.call(command, shell=False, stderr=LOG)
                LOG.write("return code = %d\n"%return_code)
                trimReads1 = os.path.basename(reads[0])
                trimReads1 = trimReads1.replace(".fastq", "_val_1.fq")
                trimReads1 = os.path.join(WorkDir, trimReads1)
                trimReads2 = os.path.basename(reads[1])
                trimReads2 = trimReads2.replace(".fastq", "_val_2.fq")
                trimReads2 = os.path.join(WorkDir, trimReads2)
                if os.path.exists(trimReads1) and os.path.exists(trimReads2):
                    comment = "trim_galore, input %s, output %s"%(reads, trimReads1+":"+trimReads2)
                    LOG.write(comment+"\n")
                    details["transformation"].append(comment)
                    platformReads[index] = trimReads1+":"+trimReads2
                else:
                    LOG.write("Problem during trim_galore: expected files not found: "+trimReads1 +", "+trimReads2)
            else:
                command.append(reads)
                LOG.write("command: "+" ".join(command))
                return_code = subprocess.call(command, shell=False, stderr=LOG)
                LOG.write("return code = %d\n"%return_code)
                trimReads = os.path.basename(reads)
                trimReads = trimReads.replace(".fastq", "_trimmed.fq")
                trimReads = os.path.join(WorkDir, trimReads)
                if os.path.exists(trimReads):
                    comment = "trim_galore, input %s, output %s"%(reads, trimReads)
                    LOG.write(comment+"\n")
                    details["transformation"].append(comment)
                    platformReads[index] = trimReads
                else:
                    raise Exception("Problem: expected file not found: "+trimReads)
    LOG.write("trim_galore trimming completed, duration = %d seconds\n"%(time()-startTrimTime))

def verifyReadPairing(readPair, output_dir):
    """ Read both files, write 3 new files: paired1, paired2, unpaired
    return paired1[:%]paired2, unpaired
    """
    LOG.write("verifyReadPairing(%s, %s)"%(readPair, output_dir))
    LOG.write("time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
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
    srf_time = time()
    LOG.write("studyReadFile() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(srf_time)), srf_time-Start_time))
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
    Read_file_type[filename] = determineReadType(lines[0], Avg_read_length[filename])
    LOG.write("found read type %s average read length %.1f\n"%(Read_file_type[filename], Avg_read_length[filename]))
    if details:
        if not "inferred_read_type" in details:
            details["inferred_read_type"] = {}
        details["inferred_read_type"][filename] = Read_file_type[filename]
        if not "estimated_read_length" in details:
            details["estimated_read_length"] = {}
        details["estimated_read_length"][filename] = Avg_read_length[filename]
    LOG.write("duration of studyReadFile was %d seconds\n"%(time() - srf_time))
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
    LOG.write("categorize_anonymous_read_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
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
                            LOG.write("Discordant fileTypes for %s vs %s, %s vs %s\n"%(filename1, filename2, Read_file_type[filename1], Read_file_type[filename2]))
                        if Read_file_type[filename1] == 'illumina':
                            if not args.illumina:
                                args.illumina = []
                            args.illumina.append(":".join(pairedFiles))
                            comment = "interpreting pair %s %s as type 'illumina'"%pairedFiles
                            LOG.write(comment+"\n")
                            details["transformation"].append(comment)
                        elif Read_file_type[filename1] == 'iontorrent':
                            if not args.iontorrent:
                                args.iontorrent = []
                            args.iontorrent.append(":".join(pairedFiles))
                            comment = "interpreting pair %s %s as type 'iontorrent'"%pairedFiles
                            LOG.write(comment+"\n")
                            details["transformation"].append(comment)
                        else: # neither illumina vs iontorrent, perhaps 'sra'
                            if Avg_read_length[filename1] < Max_short_read_length:
                                #call it illumina
                                if not args.illumina:
                                    args.illumina=[]
                                args.illumina.append(":".join(pairedFiles))
                                comment = "interpreting pair %s %s as type 'illumina'"%pairedFiles
                                LOG.write(comment+"\n")
                                details["transformation"].append(comment)
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
                comment = "interpreting file %s as type 'illumina'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
            elif Read_file_type[filename] == 'iontorrent':
                if not args.iontorrent:
                    args.iontorrent = []
                args.iontorrent.append(filename)
                comment = "interpreting file %s as type 'iontorrent'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
            elif Read_file_type[filename] == 'pacbio':
                if not args.pacbio:
                    args.pacbio = []
                args.pacbio.append(filename)
                comment = "interpreting file %s as type 'pacbio'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
            elif Read_file_type[filename] == 'nanopore':
                if not args.nanopore:
                    args.nanopore = []
                args.nanopore.append(filename)
                comment = "interpreting file %s as type 'nanopore'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
            elif Avg_read_length[filename] < Max_short_read_length:
                #call it illumina
                if not args.illumina:
                    args.illumina=[]
                args.illumina.append(filename)
                comment = "interpreting file %s as type 'illumina'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
            else:
                if not args.pacbio:
                    args.pacbio=[]
                args.pacbio.append(filename)
                comment = "interpreting file %s as type 'pacbio'"%filename
                LOG.write(comment+"\n")
                details["transformation"].append(comment)
    return

def get_runinfo(run_accession):
    """ take sra run accession (like SRR123456)
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    Altered from versionin sra_tools to handle case of multiple sra runs returned by query.
    """
    LOG.write("get_runinfo(%s)\n"%run_accession)
    runinfo = None
    runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="+run_accession
    text = urllib2.urlopen(runinfo_url).read()
    if text.startswith("Run"):
        lines = text.split("\n")
        keys   = lines[0].split(",")
        for line in lines[1:]:
            if line.startswith(run_accession):
                values = line.split(",")
                runinfo = dict(zip(keys, values))
        if not runinfo or runinfo['Platform'] not in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
            # rescue case of mis-alignment between keys and values
            LOG.write("problem in get_runinfo: sra.cgi returned:\n"+text+"\n")
            runinfo = None
            if values:
                for val in values:
                    if val in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
                        runinfo = {'Platform': val}
                        for val2 in values:
                            if val2 in ('PAIRED', 'SINGLE'):
                                runinfo['LibraryLayout'] = val2
                        break

    if not runinfo:
        LOG.write("Problem, normal runinfo request failed. Trying alternative from web page.\n")
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

def fetch_sra_files(args, details):
    """ 
    Use ftp to get all SRA files.
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    Append to appropriate parts of args (e.g., args.illumina or args.iontorrent).
    """
    LOG.write("fetch_sra_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("args.sra="+" ".join(args.sra)+"\n")
    if not args.sra:
        return
    for sra in args.sra:
        #sra_dir = ""
        #if sra.endswith(".sra"):
        #    sra_dir, sra = os.path.split(sra)
        #    sra = sra[:-4] # trim off trailing ".sra", will later detect presence of "xxx.sra" file if it exists
        sraFull = sra
        if sra.endswith(".sra"):
            sra= sra[:-4]
        runinfo = get_runinfo(sra)
        if not runinfo:
            LOG.write("runinfo for %s was empty, skipping\n"%sra)
            continue 
        LOG.write("Runinfo for %s reports platform = %s\n"%(sra, runinfo["Platform"]))

        programToUse = "fasterq-dump" # but not appropriate for pacbio or nanopore
        if runinfo['Platform'] == 'PACBIO_SMRT':
            programToUse = 'fastq-dump'
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
            programToUse = 'fastq-dump'
            if not args.nanopore:
                args.nanopore = []
            listToAddTo = args.nanopore
        else:
            LOG.write("Problem: Cannot process sra data from platform %s\n"%runinfo['Platform'])
            continue

        command = [programToUse, "-O", WorkDir, "--split-files", sraFull]
        stime = time()
        LOG.write("command = "+" ".join(command)+"\n")
        return_code = subprocess.call(command, shell=False, stderr=LOG)
        LOG.write("return_code = %d, time=%d seconds\n"%(return_code, time()-stime))
        if return_code != 0:
            LOG.write("Problem, fasterq-dump return code was %d\n"%return_code)
            LOG.write("Try one more time.\n")
            return_code = subprocess.call(command, shell=False, stderr=LOG)
            if return_code != 0:
                LOG.write("Return code on second try was %d\n"%return_code)
                LOG.write("Giving up on %s\n"%sra)
                continue

        LOG.write("Runinfo for %s reports LibraryLayout = %s\n"%(sra, runinfo['LibraryLayout']))
        fastqFiles = glob.glob(os.path.join(WorkDir, sra+"*fastq"))
        LOG.write("Number of fastq files from sra is %d\n"%len(fastqFiles))
        if len(fastqFiles):
            LOG.write("Fastq files from sra: %s\n"%(" ".join(fastqFiles)))
            if runinfo['LibraryLayout'].startswith("PAIRED") and len(fastqFiles) >= 2:
                filePair = ":".join(sorted(fastqFiles))
                listToAddTo.append(filePair)
            else:
                listToAddTo.extend(fastqFiles)
    return

def organize_read_files(args, details):
    LOG.write("\norganize_read_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
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
    if not os.path.isdir(WorkDir):
        os.mkdir(WorkDir)
    outfileName = os.path.join(WorkDir, "spades_yaml_file.txt")
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

def runQuast(contigsFile, args):
    LOG.write("runQuast() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    quastDir = os.path.join(WorkDir, "quast_out")
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

def filterContigsByMinLength(inputFile, outputFile, args, details, shortReadDepth=None, longReadDepth=None):
    """ 
    Write only sequences at or above min_length to output file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("filterContigsByMinLength(%s, %s, %d) elapsed seconds = %f\n"%(inputFile, outputFile, args.min_contig_length, time()-Start_time))
    with open(inputFile) as IN:
        with open(outputFile, 'w') as OUT:
            shortCovTags = [['sdepth', 'snorm_depth'], ['depth', 'norm_depth']][longReadDepth == None]
            longCovTags = [['ldepth', 'lnorm_depth'], ['depth', 'norm_depth']][shortReadDepth == None]
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
    details["transformation"].append(comment)

def runBandage(args, details):
    gfaFile = os.path.join(SaveDir, args.prefix+"assembly_graph.gfa")
    imageFormat = ".jpg"
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write("Bandage command =\n"+" ".join(command)+"\n")
        return_code = subprocess.call(command, shell=False, stderr=LOG)
        LOG.write("return code = %d\n"%return_code)
        if return_code == 0:
            details["Bandage plot"] = prefix+"assembly_graph.plot"+imageFormat
        else:
            LOG.write("Error creating Bandage plot\n")


def runUnicycler(args, details):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("runUnicycler: elapsed seconds = %f\n"%(time()-Start_time))
    command = ["unicycler", "-t", str(args.threads), "-o", WorkDir]
    #, "--pilon_path", args.pilon_jar]
    if args.min_contig_length:
        command.extend(("--min_fasta_length", str(args.min_contig_length)))
    command.extend(("--keep", "0")) # keep only assembly.gfa, assembly.fasta and unicycler.log
    #if args.no_pilon:
    command.append("--no_pilon")
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

    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("return code = %d\n"%return_code)

    unicyclerEndTime = time()
    elapsedTime = unicyclerEndTime - unicyclerStartTime

    details.update( { 'start_time': unicyclerStartTime,
                'end_time': unicyclerEndTime,
                'elapsed_time' : elapsedTime,
                'assembler': 'unicycler',
                'command_line': command,
                'output_path': os.path.abspath(WorkDir)
                })
    LOG.write("Duration of Unicycler run was %.1f minutes\n"%(elapsedTime/60.0))
    unicyclerLogFile = args.prefix+"unicycler.log"
    shutil.move(os.path.join(WorkDir, "unicycler.log"), os.path.join(SaveDir, unicyclerLogFile))

    if not os.path.exists(os.path.join(WorkDir, "assembly.fasta")):
        LOG.write("unicycler failed to generate assembly file.\n")
        return None
    details["contigCircular"] = []
    with open(os.path.join(WorkDir, "assembly.fasta")) as F:
        for line in F:
            if line.startswith(">"):
                details["contigCircular"].append("circular=true" in line)

    assemblyGraphFile = args.prefix+"assembly_graph.gfa"
    shutil.move(os.path.join(WorkDir, "assembly.gfa"), os.path.join(SaveDir, assemblyGraphFile))

    contigsFile = os.path.join(WorkDir, "contigs.fasta")
    shutil.move(os.path.join(WorkDir, "assembly.fasta"), contigsFile) #rename to canonical name
    return contigsFile

def runSpades(args, details):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-Start_time))
    LOG.write("runSpades: elapsed seconds = %f\n"%(time()-Start_time))

    if args.illumina and args.iontorrent:
        details["Fatal_error"] = "SPAdes cannot process both Illumina and IonTorrent reads in the same run"
        return
    command = ["spades.py", "--threads", str(args.threads), "-o", WorkDir]
    if args.singlecell:
        command.append("--sc")
    if args.iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(args)
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
        if len(single_end_reads) > 0 or len(paired_end_reads[0]) > 1:
            sys.stderr.write("SPAdes in metagenomics mode can only process a single paired-end read file\n")
            sys.exit(1);
    if args.recipe == "plasmid-spades":
        command.append("--plasmid")
    LOG.write("SPAdes command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    spadesStartTime = time()

    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("return code = %d\n"%return_code)

    spadesEndTime = time()
    elapsedTime = spadesEndTime - spadesStartTime

    details.update( { 'spades_start_time': spadesStartTime,
                'spades_end_time': spadesEndTime,
                'spades_elapsed_time' : elapsedTime,
                'assembler': 'spades',
                'command_line': command,
                'output_path': os.path.abspath(WorkDir)
                })
    elapsedTime = time() - spadesStartTime
    LOG.write("Duration of SPAdes run was %.1f hours\n"%(elapsedTime/3600.0))
    if not os.path.exists(os.path.join(WorkDir, "contigs.fasta")):
        LOG.write("spades failed to generate contigs.fasta")
        return None
    try:
        spadesLogFile = args.prefix+"spades.log"
        shutil.move(os.path.join(WorkDir, "spades.log"), os.path.join(SaveDir, spadesLogFile))
        assemblyGraphFile = args.prefix+"assembly_graph.gfa"
        shutil.move(os.path.join(WorkDir, "assembly_graph_with_scaffolds.gfa"), os.path.join(SaveDir, assemblyGraphFile))
    except Exception as e:
        LOG.write(str(e))
    contigsFile = os.path.join(WorkDir, "contigs.fasta")
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
    details["transformation"].append(comment)
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
        fastqBase = os.path.basename(read1)
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
    bamFile = runBowtie(contigFile, shortReadFastq, args.threads)
    if not bamFile:
        return None
    command = ['java', '-Xmx32G', '-jar', args.pilon_jar, '--genome', contigFile]
    if ':' in shortReadFastq:
        command.extend(('--frags', bamFile))
    else:
        command.extend(('--unpaired', bamFile))
    pilonPrefix = os.path.basename(contigFile).replace(".fasta", ".pilon")
    command.extend(('--outdir', WorkDir, '--output', pilonPrefix, '--changes'))
    command.extend(('--threads', str(args.threads)))
    tempTime = time()
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("pilon return code = %d\n"%return_code)
    LOG.write("pilon duration = %d\n"%(time() - tempTime))
    if return_code != 0:
        return None
    pilonContigs = os.path.join(WorkDir, pilonPrefix)+".fasta"
    numChanges = 0
    with open(pilonContigs.replace(".fasta", ".changes")) as CHANGES:
        for line in CHANGES:
            numChanges += 1
    LOG.write("Number of changes made by pilon was %d\n"%numChanges)

    comment = "pilon, input %s, output %s"%(contigFile, pilonContigs)
    LOG.write(comment+"\n")
    details["transformation"].append(comment)
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
    fastqBase = os.path.basename(longReadFastq)
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
    command = ["canu", "-d", WorkDir, "-p", "canu", "useGrid=false", "genomeSize=%s"%args.genome_size]
    command.extend(["maxMemory=" + str(args.memory), "maxThreads=" + str(args.threads)])
    command.append("stopOnReadQuality=false")
    """
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    if args.pacbio:
        command.append("-pacbio-raw")
        command.extend(args.pacbio) #allow multiple files
    if args.nanopore:
        command.append("-nanopore-raw")
        command.extend(args.nanopore) #allow multiple files
    LOG.write("canu command =\n"+" ".join(command)+"\n")
    #LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")

    canuStartTime = time()
    #canuLogFile = open(os.path.join(WorkDir, "canu.log"), "w")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("return code = %d\n"%return_code)
    canuEndTime = time()
    elapsedTime = canuEndTime - canuStartTime

    details.update( { 'canu_start_time': canuStartTime,
                'canu_end_time': canuEndTime,
                'canu_elapsed_time' : elapsedTime,
                'assembler': 'canu',
                'command_line': command,
                'output_path': os.path.abspath(WorkDir)
                })

    LOG.write("Duration of canu run was %.1f minutes\n"%((time() - canuStartTime)/60.0))
    contigFile = os.path.join(WorkDir, "canu.contigs.fasta")
    if not os.path.exists(contigFile):
        LOG.write("canu failed to generate contigs file.\n")
        return None
    # rename to canonical contigs.fasta
    shutil.move(os.path.join(WorkDir, "canu.contigs.fasta"), os.path.join(WorkDir, "contigs.fasta"))
    shutil.move(os.path.join(WorkDir, "canu.contigs.gfa"), os.path.join(SaveDir, prefix+"assembly_graph.gfa"))
    return os.path.join(WorkDir, "contigs.fasta")

"""
def write_html_report():
    htmlFile = os.path.join(SaveDir, args.prefix+"assembly_report.html")
    LOG.write("writing html report to %s\n"%htmlFile)
    HTML = open(htmlFile, 'w')
    HTML.write("<h1>Genome Assembly Report</h1>\n")

    HTML.write("Input reads:<br>\n")
    HTML.write("<table><tr><th>File</th><th>Platform</th><th>Layout</th<th>Avg Len</th><th>Reads</th></tr>\n")
    for item in args.illumina:
        layout = "single"
        if ":" in item:
            layout = "paired"
        HTML.write("<tr><td>"+item+"</td><td>Illumina</td><td>"+layout+"</td>")
        if item in details["read_file"
"""     

def main():
    global Start_time
    Start_time = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default="p3_assembly")
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or "%%" between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', nargs='*', help='list of IonTorrent[.gz] files or pairs, ":" between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', nargs='*', help='list of Pacific Biosciences fastq[.gz] or bam files', required=False, default=[])
    parser.add_argument('--nanopore', nargs='*', help='list of Oxford Nanotech fastq[.gz] or bam files', required=False, default=[])
    parser.add_argument('--sra', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', nargs='*', help="unspecified read files, types automatically inferred.")
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are paire-interleaved')
    parser.add_argument('--recipe', choices=['unicycler', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'auto'], help='assembler to use', default='auto')

    parser.add_argument('--racon_iterations', type=int, default=2, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=2, help='number of times to run pilon per short-read file', required=False)
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--prefix', default="", help='prefix for output files', required=False)
    parser.add_argument('--genome_size', default=Default_genome_size, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', default=1000, help='save contigs of this length or longer', required=False)
    #parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cpus', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=250, help='RAM limit in Gb')
    parser.add_argument('--trim_galore', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--pilon_jar', help='path to pilon executable or jar')
    parser.add_argument('--bandage', action='store_true', help='generate image of assembly path using Bandage')
    parser.add_argument('--params_json', help="JSON file with additional information.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    baseName = "p3x-assembly" #os.path.basename(sys.argv[0]).replace(".py", "")
    global WorkDir
    global SaveDir
    WorkDir = mkdtemp(dir=".", prefix=baseName+"_"+args.recipe+"_", suffix="_work")
    if args.outputDirectory:
        SaveDir = os.path.basename(args.outputDirectory)
    else:
        SaveDir = baseName+"_save"
    if len(args.prefix) > 0 and not args.prefix.endswith("_"):
        args.prefix += "_"
    if not os.path.exists(WorkDir):
        os.mkdir(WorkDir)
    if os.path.exists(SaveDir):
        shutil.rmtree(SaveDir)
    os.mkdir(SaveDir)
    logfileName = os.path.join(SaveDir, args.prefix + baseName + ".log")

    if args.params_json:
        parseJsonParameters(args)
    global LOG 
    sys.stderr.write("logging to "+logfileName+"\n")
    LOG = open(logfileName, 'w', 0) #unbuffered 
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(Start_time))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Temporary directory is "+WorkDir+"\n\n")
    LOG.write("Final output will be saved to "+SaveDir+"\n\n")
    details = { 'logfile' : logfileName }
    details["transformation"] = []
    if args.anonymous_reads:
        categorize_anonymous_read_files(args, details)

    if args.sra:
        fetch_sra_files(args, details)

    organize_read_files(args, details)

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
        if len(args.illumina + args.iontorrent) == 0:
            LOG.write("spades called without any short reads.\n")
        contigs = runSpades(args, details)
    elif args.recipe == "unicycler":
        if len(args.illumina + args.iontorrent + args.pacbio + args.nanopore) > 0:
            contigs = runUnicycler(args, details)
        else:
            LOG.write("unicycler called without any reads.\n")
    elif args.recipe == "canu":
        if len(args.pacbio + args.nanopore) > 0:
            contigs = runCanu(args, details)
        else:
            LOG.write("canu called without any long reads.\n")
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
                    details["transformation"].append(comment)
                else:
                    break
        
    if contigs and os.path.getsize(contigs):
        # now run pilon with each short-read file
        for i in range(0, args.pilon_iterations):
            for shortReadFastq in args.illumina + args.iontorrent:
                LOG.write("runPilon(%s, %s, args, details)\n"%(contigs, shortReadFastq))
                pilonContigFile = runPilon(contigs, shortReadFastq, args, details)
                if pilonContigFile is not None:
                    comment = "pilon: input %s, output %s"%(contigs, pilonContigFile)
                    if os.path.exists(pilonContigFile.replace(".fasta", ".changes")):
                        with open(pilonContigFile.replace(".fasta", ".changes")) as F:
                            lines = F.readlines()
                            comment += ", num_changes = %d"%len(lines)
                    LOG.write(comment+"\n")
                    details["transformation"].append(comment)
                    contigs = pilonContigFile
                else:
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
        runQuast(saveContigsFile, args)
        runBandage(args, details)

    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - Start_time)/3600))
    LOG.close()
    shutil.move(logfileName, os.path.join(args.outputDirectory, os.path.basename(logfileName)))
    fp = file(os.path.join(SaveDir, baseName+".run_details"), "w")
    json.dump(details, fp, indent=2)
    fp.close()


if __name__ == "__main__":
    main()
