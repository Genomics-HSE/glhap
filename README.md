# mt-dna
DOCUMENTATION
-------------

This software provedis you a way to calculate haplogroups depending of genotype likelihood

Dependencies
------------
C:

- htslib

https://github.com/samtools/htslib


Python

- anytree             2.8.0
- numpy               1.17.4
- pysam               0.18.0


Also you may need
PhyloTree.org-parser for making haplogrep tree in json format
https://github.com/alexeyshmelev/PhyloTree.org-parser/tree/main

How to run:
-----------
### Calculate haplogroups with samtools likelihood model:

glhap.py

you can simply use

./glhap.py array/array.json {vcf_file}.vcf

array.json is provided phylotree17 in json.

The program will produce top 10 haplogroups with their log-likelihood

### Calculate haplogroups with hamming distance model:

HOW TO USE 
g++ hamgrep.c -lhts -o hamgrep
./hamgrep f1.fa
- f1.fa is  fasta file for  called dna
YOU MUST HAVE A "filelist.txt" FILE IN FOLDER WITH EXECUTIVE FILE WHERE LIST OF FASTA FILES MUST BE WRITTEN MUST BE WRITTEN 

### Calculate haplogroups with gatk likelihood model:

g++ hamgrep.c -lhts -o hamgrep

./hamgrep {piliup}.pileup

{pileup}.pileup is pileup file, you can get it with

samtools mpileup -f ref.fa -B -o pileup.pileup  in.bam

you must have A "filelist.txt" file in folder with executive file where list of fasta files must be written. 
