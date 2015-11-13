#!/bin/bash
#################################
# 00_full_splice_site_analysis_pipeline.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-12
# ------------------
# Runs the entire splice-site analysis pipeline, from start to finish.
# Could have done this as a makefile, but I'm more familiar with Bash scripts.
# As some of the scripts can take quite a long time to run, this is probably
# best run a least over a weekend.
#
# Project folder should be set up exactly as cloned from Github.
#
# This script should be run from the root of the project folder, i.e.
# splice_sites, not splice_sites/scripts
#
#################################

############## 01_make_hg_19_fasta.sh #################
# 1) check prerequisite data
# if no hg19.fa, check for the individual files.

if [ ! -e hg19/hg19.fa ]
then
    echo "No reference genome found"
    echo "Should now run 01_make_hg19_fasta.sh. But doesn't at the moment"
fi

############## 02_get_refseq_sites.sh ##############
# Shouldn't really need to re-run this very often
# extract the RefSeq gene list and split into exons. Now removes any NR_ 
# transcripts to avoid miRNAs.
#
# output: refseq/combined.allexons - list of splice-sites ready for database
#                                    upload
# This is quite a quick script, so no need to check if it needs doing, just 
# run it?

# change to refseq folder
cd refseq_sites/

# run script
echo "DEBUG: running ../scripts/02_get_refseq_sites.sh"
#./../scripts/02_get_refseq_sites.sh

# go back to root folder
cd ..

############## 03_generate_random_sites.sh ##############

# This script takes a long time - if the files are small enough, they should be
# included in the Github repo. If not, it will need to be rerun anyway. Check
# that the output files don't exist before running. If they do, assume they are 
# up to date - only recreate if the have been deleted.

# NOTE: this is a random process, so the sites produced will not be the same 
#       each time. Unless it is possible to seed the random sort in the script?
#       NOTE: Apparently -R isn't *random*, but sorting by hash, so it should be
#             consistent, unless new data is added.

# change to random_sites folder
cd random_sites/

# check if either file doesn't exist - second check within the script will run
# for either or both files as appropriate
if [ ! -e combined.nonsplicesites ]
then
    echo "DEBUG: running ../scripts/03_generate_random_sites.sh"
    #./../scripts/03_generate_random_sites.sh
fi

# go back to root folder
cd ..

############## 04_generate_sql.sh ##############

# Now that all the data have been generated, add the variant data to a
# database.

# No later steps add any extra data to the database, so it's safe to regenerate
# from here. If anything further is to be recorded, create a new database and
# treat this one as a raw data source.

# run from within database folder, as all paths in the script are relative to it
cd database/

# NOTE: SQL generation is not yet written for random sites! Produce a few and 
#       set it up as needed.

echo "DEBUG: running ../scripts/04_generate_sql.sh"
# ./../scripts/04_generate_sql.sh

# go back to root folder
cd ..

############## 05_run_refseq_analysis.sh ##############

#TODO

# move to ther refseq fold. This script runs on both refseq and random sites,
# but handles the changing itself.
cd refseq_sites

echo "ERROR: This doesn't do anything yet."

# go back to root folder
cd ..
