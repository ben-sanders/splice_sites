#!/bin/bash

# concatenates the individual chromsome files into a whole-genome fasta

zcat chrM.fa.gz chr1.fa.gz chr2.fa.gz chr3.fa.gz chr4.fa.gz chr5.fa.gz chr6.fa.gz chr7.fa.gz chr8.fa.gz \
    chr9.fa.gz chr10.fa.gz chr11.fa.gz chr12.fa.gz chr13.fa.gz chr14.fa.gz chr15.fa.gz chr16.fa.gz \
    chr17.fa.gz chr18.fa.gz chr19.fa.gz chr20.fa.gz chr21.fa.gz chr22.fa.gz chrX.fa.gz chrY.fa.gz \
    > hg19.fa
    
# bgzip before indexing, zcat output doesn't seem to be compressed.

bgzip hg19.fa

samtools faidx hg19.fa.gz
