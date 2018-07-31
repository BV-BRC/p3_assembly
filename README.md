P3_Assembly
===========

A program to run prokaryote genome assembly for PATRIC.

Given a list of raw read files the code will run either SPAdes or Canu.

SPAdes is run if the reads are Illumina or IonTorrent (short reads).

Canu is run if the reads are PacBio or Oxford Nanopore (long reads).

If a mix of long reads and short reads is given, SPAdes will be run with long reads used to order contigs into scaffolds.





