import sys
import subprocess
import argparse
import gzip
import os.path
import re

# This script organizes a command line for SPAdes
# It auto-detects read types (illumina, iontorrent, pacbio, nanopore, fasta)
# TODO: properly handle different kinds of pacbio reads
# TODO: verify that read type identification works in general
# TODO: add read trimming option
# NOTE: I don't believe that error-correction should be be applied prior to SPAdes (which does its own version of error correction)

debug = True

def determineReadFileType(readId):
# this function should analyze sample of text from read file and return one of illumina, iontorrent, pacbio, oxfordnanopore, fasta, ...
# going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
# these patterns need to be refined and tested, some are prefixes of others, assuming numbers in the first one is based on a few observations
    if readId.startswith(">"):
        return "fasta"
    if re.match("@M\d+:\d+:\S+:\d+:\d+:\d+:\d+ \S+:\S+:\S+:\S+$", readId):
        return "illumina" # newer illumina
    if re.match("@\S+:\S+:\S+:\S+:\S+#\S+/\S+$", readId):
        return "illumina" # older illumina
    if re.match("@\S+:\S+:\S+$", readId):
        return "iontorrent" # 
# NOTE: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
    if re.match("@\S+/\S+/\S+_\S+$", readId): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
        return "pacbio_ccs_subread" # 
    if re.match("@\S+/\S+$", readId): #@<MovieName>/<ZMW_number> 
        return "pacbio_ccs" # 
    if re.match("@[a-z0-9-]+\s+read=\d+\s+ch=", readId): #based on one example, need to test more 
        return "nanopore" # 
    return "na"

def testReadIdentifiersMatch(readIds1, readIds2):
# test whether read identifiers match and are in order between to read files, presumably fastq
    if debug:
        sys.stderr.write("testReadIdentifiersMatch: \n%s \n\n%s\n\n"%(",".join(readIds1[0]), ",".join(readIds2[0])))
    limit = min(len(readIds1), len(readIds2))
    for i in range(0, limit):
        if readIds1[i] != readIds2[i]:
            if debug:
                sys.stderr.write("mismatching read[%d] ids: %s vs %s\n"%(i, readIds1[i], readIds2[i]))
            return False
    if debug:
        sys.stderr.write("all %d read ids match\n"%limit)
    return True

def studyReadFileSample(text):
    readFileType = 'na'
    readIdSample = []
    lines = text.split("\n")
    if len(lines) < 2:
        raise Exception("text sample (length %d) lacks at least 2 lines"%len(text))
    readFileType = determineReadFileType(lines[0])
    if lines[0].startswith("@"): # fastq, grab every 4th line
        for i, line in enumerate(lines):
            if i % 4 == 0:
                readIdSample.append(line.split(' ')[0]) # get part up to first space, if any 
    elif lines[0].startswith(">"): #fasta, grab every line starting with '>'
        for line in lines:
            if line.startswith(">"):
                readIdSample.append(line[1:].split(' ')[0]) # get part after '>' and up to first space, if any 
    if debug:
        sys.stderr.write("studyReadFileSample found type %s from %s\n"%(readFileType, lines[0]))
    return(readFileType, readIdSample)

def organizeReadFiles(readFileList):
    if debug:
        sys.stderr.write("organizeReadFiles: %s\n"%("\t".join(readFileList)))
    readIdSample = {}
    readFileDict = {}
    for f in readFileList:
        if f.endswith("gz"):
            F = gzip.open(f)
        else:
            F = open(f)
        text = F.read(20000) #read X number of bytes for text sample
        F.close()
        if debug:
            sys.stderr.write("  file %s:\n%s\n\n"%(f, text[0:50]))
        readIdSample[f] = []
        fileType , readIdSample[f] = studyReadFileSample(text)
        filePrefix = os.path.basename(f)
        m = re.search("(.*)_R[12]_", filePrefix) 
        if m:
            filePrefix = m.group(1)
        if fileType not in readFileDict:
            readFileDict[fileType] = {}
        if filePrefix not in readFileDict[fileType]:
            readFileDict[fileType][filePrefix] = []
        readFileDict[fileType][filePrefix].append(f)
    for t in readFileDict:
        for p in readFileDict[t]:
            if len(readFileDict[t][p]) == 2:
                testReadIdentifiersMatch(readIdSample[readFileDict[t][p][0]], readIdSample[readFileDict[t][p][1]])
    return readFileDict

def add_read_files_to_command(readFileList, command, args):
    readFileDict = organizeReadFiles(readFileList)
    if ("pacbio" in readFileDict) and (not "illumina" in readFileDict) and (not "iontorrent" in readFileDict):
        raise Exception("SPAdes cannot process PacBio reads without also either Illumina or IonTorrent reads")
    if ("illumina" in readFileDict) and ("iontorrent" in readFileDict):
        raise Exception("SPAdes cannot process both Illumina and IonTorrent reads in the same run")
    if "iontorrent" in readFileDict and not args.iontorrent:
        sys.stderr.write("Read type recognized as IonTorrent but flag not passed on command line.")
    if args.iontorrent and not "iontorrent" in readFileDict:
        sys.stderr.write("IonTorrent specified on command line but no reads of that type recognized.")
    if args.iontorrent or "iontorrent" in readFileDict:
        command += " --iontorrent" # tell SPAdes that this is the read type
    for readType in readFileDict:
        if readType == "pacbio":
            for prefix in readFileDict["pacbio"]:
                for filename in readFileDict["pacbio"][prefix]:
                    command += " --pacbio "+filname
        if readType == "illumina" or readType == "iontorrent":
            pairedEndIndex = 1
            singeReadIndex = 1
            if len(readFileDict[readType]) > 9:
                raise Exception("more than 9 input read files (or pairs), need to use Yaml file to specify that many")
            for prefix in readFileDict[readType]:
                if len(readFileDict[readType][prefix]) == 1:
                    command += " --s%d %s"%(singleReadIndex, readFileDict[readType][prefix][0])
                    singleReadIndex += 1
                elif len(readFileDict[readType][prefix]) == 2:
                    command += " --pe%d-1 %s"%(pairedEndIndex, readFileDict[readType][prefix][0])+" "+" --pe%d-2 %s"%(pairedEndIndex, readFileDict[readType][prefix][1])
                    pairedEndIndex += 1
                else:
                    raise Exception("More than two read files for prefix "+prefix)
    return command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--readFiles', nargs='+', help="whitespace sep list of read files. Types inferred automatically.")
    parser.add_argument('-o', help='output directory.', required=True)
    parser.add_argument('--debug', action = 'store_true', help='turn on debugging output', required=False)
    parser.add_argument('--trusted-contigs', help='contigs known to be good', required=False)
    parser.add_argument('--untrusted-contigs', help='contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data', required=False)
    parser.add_argument('--iontorrent', action = 'store_true', help='flag for IonTorrent data', required=False)
    parser.add_argument('--careful', action = 'store_true', help='pass careful flag to SPAdes (takes longer)', required=False)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--yamlFile', help="SPAdes 'yamlFile' specifying read files and formats. Needed for more than 9 input libraries of a single type.")
    parser.add_argument('--params', help="JSON file with additional information.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    debug = args.debug
    if debug:
        sys.stderr.write("read files: "+", ".join(args.readFiles)+"\n")
    command = "spades.py --threads %d -o %s"%(args.threads, args.o)
    if args.singlecell:
        command += " --sc"
    if args.readFiles and args.yamlFile:
        raise Exception("you must specify reads as either a read_file_list or a yamlFile, but not both")
    if args.readFiles:
        command = add_read_files_to_command(args.readFiles, command, args)
    else:
        command += " --dataset "+yamlFile
    if args.trusted_contigs:
        command += " --trusted-contigs "+args.trusted-contigs
    if args.untrusted-contigs:
        command += " --untrusted-contigs "+args.untrusted-contigs
    if args.careful:
        command += " --careful"
    if debug:
        sys.stderr.write("SPAdes command =\n"+command+"\n\n")
    if not debug:
        subprocess.call(command)


if __name__ == "__main__":
    main()
