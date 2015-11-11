#!/bin/bash
#################################
# 05_run_analysis.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Extract variant information from the database and run analysis on each 
# variant. Includes, Alamut batch, Spidex, any other tools.
# THIS VERSION IS SPECFICALLY FOR REFSEQ SITES (AND RANDOM WHEN THEY'RE DONE
#################################

# 1. Alamut batch
# ---------------

# Alamut batch wants data in tab delimited files, but not in BED format.
# From the Alamut documentation, the fields are:

# 1. Variant id (anything)
# 2. Chromosome (1-22, X, Y)
# 3. Genomic position
# 4. Reference nucleotide(s) (ACGT, or '-' for insertions)
# 5. Mutated nucleotide(s) (ACGT, or '-' for deletions)
# 6. Optional strand (1/+ or -1/-), used if --strand parameter is set to 01
# 7. Optional transcript id, used if --spectrans parameter is specified
# 8. Optional user-defined fields (e.g. heterozygosity, number of reads, etc).

# all variant data is in the database/variantdb.db, and sqlite3 database

# example query:
# "select var_id, chrom, pos_start, ref_allele, alt_allele, strand, transcript 
# from variants where source='REFSEQ'"

# then whatever output query is wanted

# run separate queries for each individual dataset? or combine as appropriate,
# e.g. REFSEQ and RANDOM?

# to select all refseq sites
# Tab separated output doesn't seem to work from the command line, so need to 
# use sed to replace | with tabs. Output is also not sorted, so sorted by chrom
# and position. This sorts by +ve chromosome order, so reverse strand genes will
# have exon order reversed. Doubt that this matters.
sqlite3 ../database/variantdb.db 'SELECT var_id, chrom, pos_start, ref_allele, alt_allele, strand, transcript FROM variants WHERE source="REFSEQ";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 \
> refseq_sites.alamut

# now the input file is formatted correctly, run alamut batch on it
# NOTE: comment out until alamut batch is actually installed!


# 2. SPIDEX
# ---------

# what is the SPIDEX input format? If it's the same as SPANR (the web-based,
# on-the fly calculated version) then it's just vcf arranged:
# CHROM POS ID  REF ALT
# with no additional data needed.
# NOTE: ID must be < 20 characters - check + correct somehow? or use db ID.
# NOTE: Check the position offset is ok, hopefully will error if ref allele 
#       doesn't match position, but it might not!

sqlite3 ../database/variantdb.db 'SELECT chrom, pos_start, var_id, ref_allele, alt_allele FROM variants WHERE source = "REFSEQ";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 \
> refseq_sites.vcf

# 3. MutPredSplice
# ----------------

# MutPredSplice is a vcf formate import, but I'm not sure if I'll actually use
# it - needs a 300+GB database setting up from scratch!

# web interface supports input in csv format:
# CHROM,POS,STRAND,REF,ALT
# or VCF as above (example file on website includes QUAL, FILTER, and INFO 
# columns, not sure if it actually NEEDS them though...)

#TODO

# 4. MESpy
# --------

# MaxEntScan is already included in Alamut, but it might be nice to try out my
# own implementation. Must be run on a per-sample basis though, so could be 
# awkward for input + output. Might need it's own wrapper script, could easily
# import + use in Python. Can account for change in score of multiple positions,
# including ranking + changes in the ranking. 

# Since this script is for refseq/random sites only, and we know where in the
# input sequence they are, restricting the input sequence to around that site 
# may help target the analysis more specifically.

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type FROM variants WHERE source="REFSEQ";' | sed s/"|"/"\t"/g | awk 'BEGIN{FS="\\t"}{print $1"\t"$2"\t"substr($3, 35, 30)"\t"substr($4, 35, 30)"\t"$5}'> refseq_sites.mespy

# this might have to be run from within the mespy folder?
python ../scripts/mespy/run_mespy.py refseq_sites.mespy

# mespy returns results as
# DB_ID VAR_ID  HIGHEST_SCORING_POSITION    RANK_CHANGE SCORE_CHANGE(%)
