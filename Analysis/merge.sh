#!/bin/bash

bam1=$1
prop1=$2
bam2=$3
prop2=$4

bam1fname=$(basename -- "$bam1")
bam1fname="${bam1fname%.*}"
bam2fname=$(basename -- "$bam2")
bam2fname="${bam2fname%.*}"

samtools view -O BAM -s $prop1 $bam1 > /tmp/bam1.bam
samtools view -O BAM -s $prop2 $bam2 > /tmp/bam2.bam
#echo $bam1 $bam2

samtools merge -o contamination_${bam1fname}_${prop1}_${bam2fname}_${prop2}.bam /tmp/bam1.bam /tmp/bam2.bam
#samtools merge -o /dev/stdout /tmp/bam1.bam /tmp/bam2.bam
rm /tmp/bam1.bam /tmp/bam2.bam
