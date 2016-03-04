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
# THIS VERSION IS SPECFICALLY FOR REFSEQ SITES (AND RANDOM WHEN THEY'RE DONE)
#################################

# 1. Alamut batch
#################

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
# database stores using +ve strand, so alamut throws an error for unexpected base (i.e complementary base) on -ve strand positions.
# awk line simply checks strand and complements -ve.

cd ../refseq_sites

sqlite3 ../database/variantdb.db 'SELECT var_id, chrom, pos_start, ref_allele, alt_allele, strand, transcript FROM variants WHERE source="REFSEQ";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 | \
awk 'BEGIN{FS="\t"; OFS="\t"}{if ($4 == "C") comp="G"; else if ($4=="G") comp="C"; else if ($4=="A") comp="T"; else comp="A"; if ($6 == "+") print $1,$2,$3,$4,$5,"1",$7; else print $1,$2,$3,comp,$5,"-1",$7}' \
> input.refseq.alamut

#TODO: run alamut (WARNING: LONG RUNNING TIME)

cd ../random_sites

sqlite3 ../database/variantdb.db 'SELECT var_id, chrom, pos_start, ref_allele, alt_allele, strand, transcript FROM variants WHERE source="RANDOM";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 | \
awk 'BEGIN{FS="\t"; OFS="\t"}{if ($4 == "C") comp="G"; else if ($4=="G") comp="C"; else if ($4=="A") comp="T"; else comp="A"; if ($6 == "+") print $1,$2,$3,$4,$5,"1",$7; else print $1,$2,$3,comp,$5,"-1",$7}' \
> input.random.alamut

#TODO: run alamut (WARNING: LONG RUNNING TIME)

cd ../literature_sites

sqlite3 ../database/variantdb.db 'SELECT var_id, chrom, pos_start, ref_allele, alt_allele, strand, transcript FROM variants WHERE source="RANDOM";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 | \
awk 'BEGIN{FS="\t"; OFS="\t"}{if ($4 == "C") comp="G"; else if ($4=="G") comp="C"; else if ($4=="A") comp="T"; else comp="A"; if ($6 == "+") print $1,$2,$3,$4,$5,"1",$7; else print $1,$2,$3,comp,$5,"-1",$7}' \
> input.literature.alamut

# 2. SPIDEX
###########

# what is the SPIDEX input format? If it's the same as SPANR (the web-based,
# on-the fly calculated version) then it's just vcf arranged:
# CHROM POS ID  REF ALT
# with no additional data needed.
# NOTE: ID must be < 20 characters - check + correct somehow? or use db ID.
# NOTE: Check the position offset is ok, hopefully will error if ref allele 
#       doesn't match position, but it might not!

cd ../refseq_sites

sqlite3 ../database/variantdb.db 'SELECT chrom, pos_start, var_id, ref_allele, alt_allele FROM variants WHERE source = "REFSEQ";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 \
> refseq.vcf

cd ../random_sites
sqlite3 ../database/variantdb.db 'SELECT chrom, pos_start, var_id, ref_allele, alt_allele FROM variants WHERE source = "RANDOM";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 \
> random.vcf

cd ../literature_sites
sqlite3 ../database/variantdb.db 'SELECT chrom, pos_start, var_id, ref_allele, alt_allele FROM variants WHERE source = "LITERATURE";' | \
sed s/'|'/'\t'/g | \
sort -k2,2 -k3n,3 \
> literature.vcf

# 3. MutPredSplice
##################

# MutPredSplice is a vcf format import, but I'm not sure if I'll actually use
# it - needs a 300+GB database setting up from scratch!

# VCF as above is fine

# 4. MESpy and MaxEntScan online
################################

# MaxEntScan is already included in Alamut, but it might be nice to try out my
# own implementation. Must be run on a per-sample basis though, so could be 
# awkward for input + output. Might need it's own wrapper script, could easily
# import + use in Python. Can account for change in score of multiple positions,
# including ranking + changes in the ranking. 

# To work out the efficiency of MES using mespy, it should only be given the
# exact splice-site sequences (or pseudo-splice) from refseq and random datasets
# this changes for 5' and 3' sites

# this cuts each type of splice site to the correct length for a single MES
# score calculation

# NOTE: giant awk command - sequence in database is from +ve strand, so any gene
#       on the negative strand needs reversing. This is now built into MESpy,
#       but needs the correct range selected from the database sequences.
#       awk checks strand, and then type, to adjust selection so the output is
#       exactly as MES wants (20I3E for 3' acceptor, 3E6I for 5' donor)

cd ../refseq_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="REFSEQ";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 48, 23)"\t"substr($4, 48, 23)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 32, 23)"\t"substr($4, 32, 23)"\t"$5"\t"$6 }}' \
> input.refseq.mespy

# split by type so can use in MaxEntScan online, and change to FASTA
cat input.refseq.mespy | cut -f2,3,6 | grep "_don" | awk 'BEGIN{FS="\t"; OFS="\n"}{print ">"$1$3,$2}' > intermediate.refseq.MESdonors.fa
python ../scripts/python/strandrevcomp.py intermediate.refseq.MESdonors.fa > input.refseq.MESdonors.fa

# run MESpy
python ../scripts/mespy/run_mespy.py input.refseq.mespy > out.refseq.mespy

# cleanup input files
mv input.refseq.mespy complete/

cd ../random_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="RANDOM";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 48, 23)"\t"substr($4, 48, 23)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 32, 23)"\t"substr($4, 32, 23)"\t"$5"\t"$6 }}' \
> input.random.mespy

# split by type so can use in MaxEntScan online, and change to FASTA
cat input.random.mespy | cut -f2,3,6 | grep "_don" | awk 'BEGIN{FS="\t"; OFS="\n"}{print ">"$1$3,$2}'  > intermediate.random.MESdonors.fa
python ../scripts/python/strandrevcomp.py intermediate.random.MESdonors.fa > input.random.MESdonors.fa

# run MESpy
python ../scripts/mespy/run_mespy.py input.random.mespy > out.random_sites.mespy

# cleanup input files
mv input.random.mespy complete/

cd ../literature_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="LITERATURE";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 48, 23)"\t"substr($4, 48, 23)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 32, 23)"\t"substr($4, 32, 23)"\t"$5"\t"$6 }}' \
> input.literature.mespy

# split by type so can use in MaxEntScan online, and change to FASTA
cat input.literature.mespy | cut -f2,3,6 | grep "_don" | awk 'BEGIN{FS="\t"; OFS="\n"}{print ">"$1$3,$2}'  > intermediate.literature.MESdonors.fa
python ../scripts/python/strandrevcomp.py intermediate.literature.MESdonors.fa > input.literature.MESdonors.fa

# run MESpy
python ../scripts/mespy/run_mespy.py input.literature.mespy > out.literature_sites.mespy

# cleanup input files
mv input.literature.mespy complete/

# mespy returns results as
# DB_ID VAR_ID  HIGHEST_SCORING_POSITION    RANK_CHANGE SCORE_CHANGE(%)
# score change is percentage relative to wildtype.
# increased score = positive percentage (e.g. 10%)
# decreased score = negative percentage
# 100% decrease would be variant score of 0. -150% would be going to negative
# half the original (e.g. 10 to -5)
# therefore according to BPGs, any variant with score change -10 or below is 
# categorised as splice affecting. Analyse this later, rather than adding an awk
# step above - can then try different cutoffs to see the effect - ROC curves
# Actually, given that, do I even want the percentage auto calculated?
# that might be best, as the other tools will give different scores (e.g. 0-100,
# 0.0 - 1.0, etc)

# 5. GeneSplicer
################

# GeneSplicer takes fast inputs, but doesn't seem to account for multiple
# samples - it concatenates them all together!?
# need to feed it individual files, which obviously would be a pain so will
# extract database as one file and then pull out the individual sequences for 
# GeneSplicer

# will run this in a separate script, just use this for the initial database
# extraction.

#TODO

# 6. PWM
########

# Check against the derived PWMs from the S&S analysis subsection
# input format is very similary to MESpy.

cd ../refseq_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="REFSEQ";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 50, 15)"\t"substr($4, 50, 15)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 38, 15)"\t"substr($4, 38, 15)"\t"$5"\t"$6 }}' \
> input.refseq.pwm

python ../scripts/python/pwm_score.py input.refseq.pwm > out.refseq.pwm

# tidy up
mv input.refseq.pwm complete/

cd ../random_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="RANDOM";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 50, 15)"\t"substr($4, 50, 15)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 38, 15)"\t"substr($4, 38, 15)"\t"$5"\t"$6 }}' \
> input.random.pwm

python ../scripts/python/pwm_score.py input.random_sites.pwm > out.random.pwm

# tidy up
mv input.random.pwm complete/

cd ../literature_sites

sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="LITERATURE";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 50, 15)"\t"substr($4, 50, 15)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 38, 15)"\t"substr($4, 38, 15)"\t"$5"\t"$6 }}' \
> input.literature.pwm

python ../scripts/python/pwm_score.py input.literature_sites.pwm > out.literature.pwm

# tidy up
mv input.literature.pwm complete/

# 7. NNSPLICE online FASTA
##########################

# For generic use with some web tools

cd ../refseq_sites

# get the sequence info and process to the right lengths
sqlite3 ../database/variantdb.db 'SELECT db_id, var_id, ref_seq, alt_seq, type, strand FROM variants WHERE source="REFSEQ";' | \
sed s/"|"/"\t"/g | \
awk 'BEGIN{FS="\\t"}{if ($6 == "-") { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 46, 9)"\t"substr($4, 46, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 48, 23)"\t"substr($4, 48, 23)"\t"$5"\t"$6 } else if ($6 == "+" ) { if ($5 == "don" ) print $1"\t"$2"\t"substr($3, 48, 9)"\t"substr($4, 48, 9)"\t"$5"\t"$6;  else print $1"\t"$2"\t"substr($3, 32, 23)"\t"substr($4, 32, 23)"\t"$5"\t"$6 }}' \
> input.refseq.nnsplice.fa

#########TODO : FINISH