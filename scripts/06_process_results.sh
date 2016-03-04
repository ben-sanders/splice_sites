#!/bin/bash
#################################
# 05_combine_results.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Integrates the output from the various analysis tools (Alamut, Spidex, etc),
# works out the changes between wt and variant, and runs scripts to compare +
# analyse the whole dataset.
#################################

################ arrange results ######################
# for each tool, just want ID, WTSCORE, and VARSCORE
# this will make analysis with R easier
# REMEMBER TO NORMALISE SCORES TO 0-1

#### RefSeq
###########
cd ../refseq_sites/

# MESPY output
cat out.refseq.mespy | \
cut -f2,5,6 \
> scores.refseq.mespy

# PWM output
cat out.refseq.pwm | \
cut -f1,4,5 \
> scores.refseq.pwm

# GeneSplicer output is lacking var score - look into how it's run!
cp out.refseq.genesplicer scores.refseq.genesplicer

# Alamut output
# SSFL
cat out.refseq.alamut | cut -f1,33,38 > scores.refseq.ssfl.alamut

# MaxEntScan
cat out.refseq.alamut | cut -f1,34,39 > scores.refseq.maxentscan.alamut

# NNSplice
cat out.refseq.alamut | cut -f1,35,40 > scores.refseq.nnsplice.alamut

# GeneSplicer
cat out.refseq.alamut | cut -f1,36,41 > scores.refseq.genesplicer.alamut

#HSF
cat out.refseq.alamut | cut -f1,37,42 > scores.refseq.ssfl.alamut

# tidy up
mv out.refseq.* outputs/
mw scores.refseq.* scores/

#### random sites
#################

# copy from above + change refseq to random.

################ R script and plotting ########################
#cd ../refseq_sites

#R < ../scripts/r/MES_plots.r --no-save
#R < ../scripts/r/genesplicer_plots.r --no-save