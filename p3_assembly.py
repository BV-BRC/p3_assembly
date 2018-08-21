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

def runTrimmomatic(args):
    if args.illumina:
        if not os.path.isfile(args.illuminaAdapters): #default illuminaAdapters is "illumina_adapters.fa", should be in same directory as this script
            scriptDir = os.path.dirname(sys.argv[0])
            if os.path.isfile(os.path.join(scriptDir, args.illuminaAdapters)):
                args.illuminaAdapters = os.path.join(scriptDir, args.illuminaAdapters)
        trimmedIllumina = []
        for item in args.illumina:
            if ':' in item or '%' in item:
            # paired-end read
                separator = ':'
                if ':' in item:
                    read1, read2 = item.split(':')
                else:
                    read1, read2 = item.split('%')
                    separator = '%'
                read1OutBase = os.path.join(args.output_dir, os.path.basename(read1))
                if read1OutBase.endswith(".fq.gz"):
                    read1OutBase = read1OutBase[:-6]
                elif read1OutBase.endswith(".fq"):
                    read1OutBase = read1OutBase[:-3]
                elif read1OutBase.endswith(".fastq"):
                    read1OutBase = read1OutBase[:-6]
                elif read1OutBase.endswith(".fastq.gz"):
                    read1OutBase = read1OutBase[:-9]
                elif read1OutBase.endswith(".gz"):
                    read1OutBase = read1OutBase[:-3]
                read1OutBase += "_trim"
                read2OutBase = os.path.join(args.output_dir, os.path.basename(read2))
                if read2OutBase.endswith(".fq.gz"):
                    read2OutBase = read2OutBase[:-6]
                elif read2OutBase.endswith(".fq"):
                    read2OutBase = read2OutBase[:-3]
                elif read2OutBase.endswith(".fastq"):
                    read2OutBase = read2OutBase[:-6]
                elif read2OutBase.endswith(".fastq.gz"):
                    read2OutBase = read2OutBase[:-9]
                elif read2OutBase.endswith(".gz"):
                    read2OutBase = read2OutBase[:-3]
                read2OutBase += "_trim"
                
                trimLog = read1OutBase+".log"
                #command = "java -jar %s PE -threads %d -trimlog %s "%(args.pathToTrimmomatic, args.threads, trimLog)
                #command += " %s %s %s %s %s %s "%(read1, read2, read1OutBase+"_P.fq", read1OutBase+"_U.fq", read2OutBase+"_P.fq", read2OutBase+"_U.fq")
                #command += " ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:%d:%d"%(args.illuminaAdapters, args.trimmomaticWindow, args.trimmomaticMinQual)
                command = ["java", "-jar", args.pathToTrimmomatic, "PE", "-threads", str(args.threads), "-trimlog", trimLog]
                command.extend([read1, read2, read1OutBase+"_P.fq", read1OutBase+"_U.fq", read2OutBase+"_P.fq", read2OutBase+"_U.fq"])
                if os.path.isfile(args.illuminaAdpaters):
                    command.append("ILLUMINACLIP:%s:2:30:10"%args.illuminaAdapters)
                command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
                if args.debug:
                    sys.stderr.write("command = "+" ".join(command)+"\n")
                subprocess.call(command, shell=False)
                trimmedIllumina.append(separator.join(read1OutBase+"_P.fq", read2OutBase+"_P.fq"))
                trimmedIllumina.append(read1OutBase+"_U.fq")
                trimmedIllumina.append(read2OutBase+"_U.fq")

            else:
                # an unpaired read file
                outBase = os.path.join(args.output_dir, os.path.basename(item))+"_trim"
                trimLog = outBase+".log"
                #command = "java -jar %s SE -threads %d -trimlog %s "%(args.pathToTrimmomatic, args.threads, trimLog)
                #command += " %s %s "%(item, outBase)
                #command += " ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:%d:%d"%(args.illuminaAdapters, args.trimmomaticWindow, args.trimmomaticMinQual)
                command = ["java", "-jar", args.pathToTrimmomatic, "SE", "-threads", str(args.threads), "-trimlog", trimLog]
                command.extend([item, outBase])
                if os.path.isfile(args.illuminaAdpaters):
                    command.append("ILLUMINACLIP:%s:2:30:10"%args.illuminaAdapters)
                command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
                if args.debug:
                    sys.stderr.write("command = "+" ".join(command)+"\n")
                subprocess.call(command, shell=False)
                trimmedIllumina.append(outBase)
        args.illumina = trimmedIllumina # replace original list of illumina reads with trimmed versions

    if args.iontorrent:
        trimmedIontorrent = []
        for item in args.iontorrent:
            if ':' in item:
            # paired-end read
                read1, read2 = item.split(':')
                read1OutBase = os.path.join(args.output_dir, os.path.basename(read1))
                if read1OutBase.endswith(".fq.gz"):
                    read1OutBase = read1OutBase[:-6]
                elif read1OutBase.endswith(".fq"):
                    read1OutBase = read1OutBase[:-3]
                elif read1OutBase.endswith(".fastq"):
                    read1OutBase = read1OutBase[:-6]
                elif read1OutBase.endswith(".fastq.gz"):
                    read1OutBase = read1OutBase[:-9]
                elif read1OutBase.endswith(".gz"):
                    read1OutBase = read1OutBase[:-3]
                read1OutBase += "_trim"
                read2OutBase = os.path.join(args.output_dir, os.path.basename(read2))
                if read2OutBase.endswith(".fq.gz"):
                    read2OutBase = read2OutBase[:-6]
                elif read2OutBase.endswith(".fq"):
                    read2OutBase = read2OutBase[:-3]
                elif read2OutBase.endswith(".fastq"):
                    read2OutBase = read2OutBase[:-6]
                elif read2OutBase.endswith(".fastq.gz"):
                    read2OutBase = read2OutBase[:-9]
                elif read2OutBase.endswith(".gz"):
                    read2OutBase = read2OutBase[:-3]
                read2OutBase += "_trim"
                
                trimLog = read1OutBase+".log"
                
                command = ["java", "-jar", args.pathToTrimmomatic, "PE", "-threads", str(args.threads), "-trimlog", trimLog]
                command.extend([read1, read2, read1OutBase+"_P.fq", read1OutBase+"_U.fq", read2OutBase+"_P.fq", read2OutBase+"_U.fq"])
                command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
                if args.debug:
                    sys.stderr.write("command = "+" ".join(command)+"\n")
                subprocess.call(command, shell=False)
                trimmedIontorrent.append(":".join(read1OutBase+"_P.fq", read2OutBase+"_P.fq"))
                trimmedIontorrent.append(read1OutBase+"_U.fq")
                trimmedIontorrent.append(read2OutBase+"_U.fq")

            else:
            # not paired
                outBase = os.path.join(args.output_dir, os.path.basename(item))+"_trim"
                trimLog = outBase+".log"
                command = ["java", "-jar", args.pathToTrimmomatic, "SE", "-threads", str(args.threads), "-trimlog", trimLog]
                command.extend([item, outBase])
                command.append("SLIDINGWINDOW:%d:%d"%(args.trimmomaticWindow, args.trimmomaticMinQual))
            if args.debug:
                sys.stderr.write("command = "+" ".join(command)+"\n")
            subprocess.call(command, shell=False)
            trimmedIontorrent.append(outBase)
            args.iontorrent = trimmedIontorrent # replace original list of iontorrent reads with trimmed versions

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
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    outfileName = os.path.join(args.output_dir, "spades_yaml_file.txt")
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
    #if ("illumina_pe" in args.output_dirr "illumina_se" in args) and ("iontorrent_pe" in args or "iontorrent_se" in args):
    if args.illumina and args.iontorrent:
        raise Exception("SPAdes cannot process both Illumina and IonTorrent reads in the same run")
    testPairedReadIdentifiersMatch(args)
    command = ["spades.py", "--threads", str(args.threads), "-o", args.output_dir]
    if args.singlecell:
        command.append("--sc")
    if args.iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(args)
    command.extend(["--dataset", yamlFile])
    if args.trusted_contigs:
        command.extend(["--trusted-contigs", args.trusted-contigs])
    if args.untrusted_contigs:
        command.extend(["--untrusted-contigs", args.untrusted-contigs])
    if args.careful:
        command.append("--careful")
    if args.debug:
        sys.stderr.write("SPAdes command =\n"+" ".join(command)+"\n")
        sys.stderr.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    if not args.debug:
        subprocess.call(command, shell=False)
    if args.quast:
        if not args.quast_path:
            args.quast_path = os.path.dirname(sys.argv[0]) # try directory this script is in
        subprocess.call([os.path.join(args.quast_path, "quast.py"), "-o", "quast_out", "--gene-finding", "contigs.fasta", "scaffolds.fasta"], shell=False)

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
    command = ["canu", "-d", args.output_dir, "-p", args.canu_prefix, "gnuplotTested=true", "useGrid=false", "genomeSize=%s"%args.genome_size]
    if args.pacbio:
        command.append("-pacbio-raw")
        command.extend(args.pacbio) #allow multiple files
    if args.nanopore:
        command.append("-nanopore-raw")
        command.extend(args.nanopore) #allow multiple files
    if args.debug:
        sys.stderr.write("canu command =\n"+" ".join(command)+"\n")
        sys.stderr.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    if not args.debug:
        subprocess.call(command, shell=False)
    if args.quast:
        if not args.quast_path:
            args.quast_path = os.path.dirname(sys.argv[0]) # try directory this script is in
        subprocess.call([os.path.join(args.quast_path, "quast.py"), "-o", "quast_out", "--gene-finding", args.canu_prefix+".contigs.fasta", args.canu_prefix+".unitigs.fasta"], shell=False)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output_dir', default='.', help='output directory.', required=False)
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or "%%" between mate-pairs', required=False)
    illumina_or_iontorrent.add_argument('--iontorrent', nargs='*', help='list of IonTorrent[.gz] files or pairs, ":" between paired-end-files', required=False)
    parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data', required=False)
    parser.add_argument('--pacbio', nargs='*', help='list of Pacific Biosciences fastq[.gz] or bam files "," between libraries', required=False)
    parser.add_argument('--nanopore', nargs='*', help='list of Oxford Nanotech fastq[.gz] or bam files "," between libraries', required=False)
    parser.add_argument('--canu_prefix', default='canu', help='prefix for canu output', required=False)
    parser.add_argument('--genome_size', default=default_genome_size, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--anonymous_reads', nargs='*', help="unspecified read files, types automatically inferred.")
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('--careful', action = 'store_true', help='pass careful flag to SPAdes (takes longer)', required=False)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--bytes_to_sample', type=int, default=default_bytes_to_sample, help='how much to sample from read files to test file type')
    parser.add_argument('--runTrimmomatic', action = 'store_true', help='run trimmomatic on Illumina or Iontorrent fastq files')
    parser.add_argument('--trimmomatic_jar', default='trimmomatic.jar', help='trimmomatic jar file, with path')
    parser.add_argument('--illuminaAdapters', default='illumina_adapters.fa', help='illumina adapters file, looked for in script directory')
    parser.add_argument('--trimmomaticWindow', type=int, default=4, help='window width for trimming')
    parser.add_argument('--trimmomaticMinQual', type=int, default=15, help='min score of window below which 3\' end is trimmed')
    parser.add_argument('--quast', action = 'store_true', help='run quast for assembly quality statistics')
    parser.add_argument('--quast_path', help='path to quast.py (excluding script name)')
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
    if args.runTrimmomatic:
        runTrimmomatic(args)
# if any illumina or iontorrent reads present, must use SPAdes (long-reads can be present), else use canu for long-reads
    if args.illumina or args.iontorrent:
        runSpades(args)
    else:
        runCanu(args)

if __name__ == "__main__":
    main()
