#!/bin/bash
counter=1
while IFS= read -r url; do
  samtools view -b "$url"  chrY > ydna/in${counter}.bam
  samtools view -b "$url"  chrM > mtdna/in${counter}.bam
done < files.txt