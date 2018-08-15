import sys
import subprocess
import argparse
import gzip
import os
import os.path
import re

# This script organizes a command line for either SPAdes or canu as appropriate (canu if pacbio or nanopore, spades if illumina or iontorrent (possibly plus pacbio or nanopore))
# It auto-detects read types (illumina, iontorrent, pacbio, nanopore, fasta)
# TODO: properly handle different kinds of pacbio reads
# TODO: verify that read type identification works in general
# TODO: add read trimming option

default_genome_size = "5m"
default_bytes_to_sample = 20000
readIdSample = {}
readFileType = {}
maxReadLength = {}

def determineReadFileType(readId):
# this function should analyze sample of text from read file and return one of illumina, iontorrent, pacbio, oxfordnanopore, fasta, ...
# going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
# these patterns need to be refined and tested, some are prefixes of others, assuming numbers in the first one is based on a few observations
    if readId.startswith(">"):
        return "fasta"
    #@D00553R:173:HG53VBCXY:2:1101:1235:2074 1:N:0:ACAGTGAT
    if re.match("@[A-Z]\S+:\d+:\S+:\d+:\d+:\d+:\d+ \S+:\S+:\S+:\S+$", readId):
        return "illumina" # newer illumina
    if re.match("@\S+:\S+:\S+:\S+:\S+#\S+/\S+$", readId):
        return "illumina" # older illumina
    if re.match("@\S+:\S+:\S+$", readId):
        return "iontorrent" # 
    if re.match("@[SED]RR\d+\.\d+", readId):
        return "sra" # 
# NOTE: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
    if re.match("@\S+/\S+/\S+_\S+$", readId): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
        return "pacbio_ccs_subread" # 
    if re.match("@\S+/\S+$", readId): #@<MovieName>/<ZMW_number> 
        return "pacbio_ccs" # 
#@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
    if re.match("@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", readId): #based on one example, need to test more 
        return "nanopore" # 
    return "na"

def studyReadFile(filename):
    readFileType[filename] = 'na'
    readIdSample[filename] = []
    if filename.endswith("gz"):
        F = gzip.open(filename)
    else:
        F = open(filename)
    text = F.read(default_bytes_to_sample) #read X number of bytes for text sample
    F.close()
    if args.debug:
        sys.stderr.write("  file %s:\n%s\n\n"%(filename, text[0:50]))
    readIdSample[filename] = []
    lines = text.split("\n")
    if len(lines) < 2:
        raise Exception("text sample (length %d) lacks at least 2 lines"%len(text))
    readFileType[filename] = determineReadFileType(lines[0])
    read_format = None
    if lines[0].startswith("@"):
        read_format = 'fastq'
    elif lines[0].startswith(">"):
        read_format = 'fasta'
    maxReadLength[filename] = 0
    if read_format == 'fastq':
        for i, line in enumerate(lines):
            if i % 4 == 0:
                readIdSample.append(line.split(' ')[0]) # get part up to first space, if any 
            elif i % 4 == 1:
                maxReadLength[filename] = max(maxReadLength[filename], len(line)-1)
    elif read_format == 'fastq':
        read = ''
        for line in lines:
            if line.startswith(">"):
                readIdSample.append(line[1:].split(' ')[0]) # get part after '>' and up to first space, if any 
                maxReadLength[filename] = max(maxReadLength[filename], len(read))
                read = ''
            else:
                read += line.rstrip()
    if args.debug:
        sys.stderr.write("studyReadFile found type %s max read length %d\nfrom %s\n"%(readFileType[filename], maxReadLength[filename], lines[0]))
    return(readFileType[filename], readIdSample[filename], maxReadLength[filename])

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

def categorize_anonymous_read_files(args):
    if args.debug:
        sys.stderr.write("categorize_anonymous_read_files: %s\n"%("\t".join(args.anonymous_reads)))
    for filename in args.anonymous_reads:
        fileType[filename], readIdSample[filename], maxReadLength = studyReadFile(filename)

    # try to find paired files
    membersOfPairs = set()
    for i, filename1 in enumerate(args.anonymous_readsi[:-1]):
        for filename2 in args.anomymousReads[i+1:]:
            if fileType[filename1] != fileType[filename2]:
                continue
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
                        if fileType[filename1] == 'illumina':
                            if not args.illumina:
                                args.illumina = []
                            args.illumina_pe.append(":".join(pairedFiles))
                        if fileType[filename1] == 'iontorrent':
                            if not args.args.iontorrent:
                                args.iontorrent = []
                            args.iontorrent.append(":".join(pairedFiles))
                        membersOfPairs.add(pairedFiles[0])
                        membersOfPairs.add(pairedFiles[1])
                except:
                    pass

    for filename in args.anonymousFiles:
        if filename not in membersOfPairs:
            if fileType[filename1] == 'illumina':
                if not args.illumina:
                    args.illumina = []
                args.illumina_se.append(filename)
            if fileType[filename1] == 'iontorrent':
                if not args.iontorrent:
                    args.iontorrent = []
                args.iontorrent.append(filename)
    return

def testPairedReadIdentifiersMatch(args):
    if args.debug:
        sys.stderr.write("testPairedReadIdentifiersMatch: \n")
    all_read_files = []
    if args.illumina:
        all_read_files = args.illumina
    elif args.iontorrent:
        all_read_files = args.iontorrent
    for item in all_read_files:
        if ":" in item:
            filePair = item.split(":")
        elif "%" in item:
            filePair = item.split("%") # illumina mate-filePair
        else:
            continue
        for f in filePair:
            if f not in readIdSample:
                studyReadFile(f)
                #if readFileType[f] not in t:
                #    raise Exception("Hey! filetype conflict for %s: looks like %s, but input as %s"%(f, fileType, t))
        for readIdPair in zip(readIdSample[filePair[0]], readIdSample[filePair[1]]):
            if readIdPair[0] != readIdPair[1]:
                raise Exception("Hey! read IDs do not match for files %s and %s"%filePair)
    return True

def writeSpadesYamlFile(args):
    if not os.path.isdir(args.o):
        os.mkdir(args.o)
    outfileName = os.path.join(args.o, "spades_yaml_file.txt")
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    single_end_reads = []
    paired_end_reads = [[], []]
    mate_pair_reads = [[], []]
    all_read_files = []
    if args.illumina:
        all_read_files = args.illumina
    elif args.iontorrent:
        all_read_files = args.iontorrent
    for item in all_read_files:
        if ":" in item:
            pair = item.split(":")
            f = os.path.abspath(pair[0])
            paired_end_reads[0].append(f)
            f = os.path.abspath(pair[1])
            paired_end_reads[1].append(f)
        elif "%" in item:
            pair = item.split("%")
            f = os.path.abspath(pair[0])
            mate_pair_reads[0].append(f)
            f = os.path.abspath(pair[1])
            mate_pair_reads[1].append(f)
        else:
            single_end_reads.append(os.path.abspath(item))
    if len(single_end_reads):
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(single_end_reads))
        OUT.write("\"\n    ]\n  }\n")
    if len(paired_end_reads[0]):
        OUT.write("  {\n    orientation: \"fr\",\n")
        OUT.write("    type: \"paired-end\",\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(paired_end_reads[0]))
        OUT.write("\"\n    ],\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(paired_end_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
    if len(mate_pair_reads[0]):
        OUT.write("  {\n    orientation: \"rf\",\n")
        OUT.write("    type: \"mate-pairs\",\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[0]))
        OUT.write("\"\n    ]\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
    if args.pacbio:
        pacbio_reads = []
        for f in args.pacbio:
            pacbio_reads.append(os.path.abspath(f))
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
    if args.nanopore:
        nanopore_reads = []
        for f in args.nanopore:
            nanopore_reads.append(os.path.abspath(f))
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")

    OUT.write("]\n")
    OUT.close()
    return(outfileName)    

def runSpades(args):
    #if ("illumina_pe" in args or "illumina_se" in args) and ("iontorrent_pe" in args or "iontorrent_se" in args):
    if args.illumina and args.iontorrent:
        raise Exception("SPAdes cannot process both Illumina and IonTorrent reads in the same run")
    testPairedReadIdentifiersMatch(args)
    command = "spades.py --threads %d -o %s"%(args.threads, args.o)
    if args.singlecell:
        command += " --sc"
    if args.iontorrent:
        command += " --iontorrent" # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(args)
    command += " --dataset " + yamlFile
    if args.trusted_contigs:
        command += " --trusted-contigs "+args.trusted-contigs
    if args.untrusted_contigs:
        command += " --untrusted-contigs "+args.untrusted-contigs
    if args.careful:
        command += " --careful"
    if args.debug:
        sys.stderr.write("SPAdes command =\n"+command+"\n")
        sys.stderr.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    if not args.debug:
        subprocess.call(command, shell=True)

def runCanu(args):
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
    tempdir = 'XXX'
    prefix = 'YYY'
    command = "canu -d %s -p %s gnuplotTested=true useGrid=false genomeSize=%s "%(tempdir, prefix, args.genome_size)
    if args.pacbio:
        command += " -pacbio-raw " + " ".join(args.pacbio) #allow multiple files
    if args.nanopore:
        command += " -nanopore-raw " + " ".join(args.nanopore) #allow multiple files
    if args.debug:
        sys.stderr.write("canu command =\n"+command+"\n")
        sys.stderr.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    if not args.debug:
        subprocess.call(command, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='output directory.', required=True)
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or "%" between mate-pairs', required=False)
    illumina_or_iontorrent.add_argument('--iontorrent', nargs='*', help='list of IonTorrent[.gz] files or pairs, ":" between paired-end-files', required=False)
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data', required=False)
    parser.add_argument('--pacbio', nargs='*', help='list of Pacific Biosciences fastq[.gz] or bam files "," between libraries', required=False)
    parser.add_argument('--nanopore', nargs='*', help='list of Oxford Nanotech fastq[.gz] or bam files "," between libraries', required=False)
    parser.add_argument('--genome_size', default=default_genome_size, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--anonymous_reads', nargs='*', help="unspecified read files, types automatically inferred.")
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('--careful', action = 'store_true', help='pass careful flag to SPAdes (takes longer)', required=False)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--bytes_to_sample', type=int, default=default_bytes_to_sample, help='how much to sample from read files to test file type')
    parser.add_argument('--debug', action = 'store_true', help='turn on debugging output', required=False)
    #parser.add_argument('--params', help="JSON file with additional information.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    if args.debug:
        print "args= "+str(args)
    if args.anonymous_reads:
        categorize_anonymous_read_files(args)
    # if any illumina or iontorrent reads present, must use SPAdes (long-reads can be present), else use canu for long-reads
    if args.illumina or args.iontorrent:
        runSpades(args)
    else:
        runCanu(args)


if __name__ == "__main__":
    main()
