# glhap
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


# gl_cont

This software provides you to estimate contamination in DNA based on quality score.

Dependencies
------------
- simlord                   1.0.4 
- numpy                     1.22.3 
- tqdm                      4.64.0 
- biopython                 1.78
- cython                    0.29.30
- pysam                     0.19.1



How to run:
-----------
for first use you should type
python setup.pu build_ext  --inplace

then

./gl_cont.py bam.bam ref.ref contaminants.fa nIter
where
nIter - number of iteration for MCMC


### Simulate data
------------------
You may generate data with script based on simlord

python contamsim.py hap1.fa prop1 [...] hapn.fa propn
prop1 is number of reads divided by 100000   
