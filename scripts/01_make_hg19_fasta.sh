#!/bin/bash
#################################
# 01_make_hg19_fasta.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# from a given RefSeq gene list (refseq_sites/gene_list/refgenes.txt - this may
# get replaced with a download link, to ensure it's always up to date), extracts
# splice sites BED file 
# The gene list location if defined in the python script, so that would need to 
# be edited if the script was also going to download it. Probably good practice 
# anyway, to avoid hardcoding too much.
#################################

# this is all assuming these files exist - how can I check a list?
# for i in *.fa.gz would work for existing files...

cd ../hg19

# concatenates the individual chromsome files into a whole-genome fasta

zcat chrM.fa.gz chr1.fa.gz chr2.fa.gz chr3.fa.gz chr4.fa.gz chr5.fa.gz \
chr6.fa.gz chr7.fa.gz chr8.fa.gz chr9.fa.gz chr10.fa.gz chr11.fa.gz \
chr12.fa.gz chr13.fa.gz chr14.fa.gz chr15.fa.gz chr16.fa.gz chr17.fa.gz \
chr18.fa.gz chr19.fa.gz chr20.fa.gz chr21.fa.gz chr22.fa.gz chrX.fa.gz \
chrY.fa.gz > hg19.fa

# index the whole-genome fasta
samtools faidx hg19.fa

# the index also works as the genome file for some later bedtools functions
# since it's small and makes it clearer, copy this as hg19.genome
cp hg19.fa.fai hg19.genome
    
# unfortunately bedtools can't read gzipped fasta, so can't compress and save
# some disk space.
