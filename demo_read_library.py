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
from ReadLibrary import ReadLibrary
MAX_BASES=1e11

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='demo_read_libary_work')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    parser.add_argument('--fasta', metavar='files', nargs='*', help='list of fasta[.gz] files', required=False, default=[])
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--max_bases', type=int, default=MAX_BASES, help='downsample reads if more than this total bases.')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=4)
    parser.add_argument('--trim', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--normalize', action='store_true', help='normalize read depth with BBNorm at default settings')

    args = parser.parse_args()
    
    WORK_DIR = args.outputDirectory
    ReadLibrary.WORK_DIR = WORK_DIR   
    if not os.path.exists(WORK_DIR):
        os.mkdir(WORK_DIR)
    ReadLibrary.NUM_THREADS = args.threads
    ReadLibrary.LOG = sys.stderr
    read_list = []
    if args.anonymous_reads:
        for item in args.anonymous_reads:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, work_dir=WORK_DIR)
            else:
                readLib = ReadLibrary(item, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.illumina:
        platform = 'illumina'
        for item in args.illumina:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, platform=platform, work_dir=WORK_DIR)
            else:
                interleaved = args.interleaved and item in args.interleaved
                readLib = ReadLibrary(item, platform=platform, interleaved=interleaved, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.iontorrent:
        platform = 'iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            readLib = ReadLibrary(item, platform=platform, interleaved=interleaved, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.pacbio:
        platform='pacbio'
        for item in args.pacbio:
            readLib = ReadLibrary(item, platform=platform, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.nanopore:
        print("process nanopore items")
        platform='nanopore'
        for item in args.nanopore:
            print(item)
            readLib = ReadLibrary(item, platform=platform, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.fasta:
        for item in args.fasta:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, work_dir=WORK_DIR)
            else:
                readLib = ReadLibrary(item, work_dir=WORK_DIR)
            read_list.append(readLib)

    # move into working directory so that all files are local
    os.chdir(WORK_DIR)

    if args.trim:
        for read_set in read_list:
            if read_set.length_class == "short" and read_set.format == 'fastq': # TrimGalore only works on short fastq reads
                read_set.trim_short_reads()

    if args.normalize:
        for read_set in read_list:
            read_set.normalize_read_depth()

    if args.max_bases:
        for read_set in read_list:
            if read_set.num_bases > args.max_bases:
                read_set.down_sample_reads(args.max_bases)

    htmlFile = os.path.join("demo_read_library.html")
    HTML = open(htmlFile, 'wt')
    HTML.write("<!DOCTYPE html><html>\n")

    HTML.write('<div>\n<h2>Preprocessing of Reads</h2>\n')
    for read_set in read_list:
        read_set.writeHtmlSection(HTML)
    HTML.write("</div>\n")

    HTML.write("</html>\n")

if __name__ == '__main__':
    main()
