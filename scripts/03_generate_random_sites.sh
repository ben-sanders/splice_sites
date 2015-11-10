#!/bin/bash
################################
# 03_generate_random_sites.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Select locations of AG/GT dinucleotides (as in splice-sites) that are know
# not to be functional splice sites. Inverts the list of known sites to get
# regions that are known not to contain splice sites, and intersects this with
# genome-wide lists of AG and GT dinucleotide positions. Then shuffles the list
# and selects a matched number of each site type to the number in the refseq
# list.
#################################

# script variables - list of splice-site positions, number of positions, and
# number of each type (i.e. total number / 2)
EXONFILE="../refseq_sites/combined.allexons"
SITENUM=`wc -l "$EXONFILE" | awk '{print $1}'`
PERTYPE=$(( $SITENUM / 2 ))

#NOTE: This takes a long time to run, as it has to iterate over the entire
#      human genome several times

# bedtools complement reverses the selection - returns the regions of the
# genome file that aren't covered by intervals in the combined.allexons bedfile.
cat "$EXONFILE" | bedtools slop -i stdin -g ../hg19/hg19.genome -b 20 | \
bedtools complement -i stdin -g ../hg19/hg19.genome > inverse.allexons

# now i need to go through these ranges, select somewhere from within them, and
# return those selections as another bed file, to get ref and generate alt 
# alleles

# ALTERNATIVE: generate bed file of all AG/GT sites in the genome, get rid of
# all that intersect with slopped combined.allexons ranges, and select at random
# ensure a matching split between donors and acceptors.

python ../scripts/nonsplicesites.py ../hg19/hg19.fa AG | \
bedtools intersect -sorted -a inverse.allexons -b stdin | sort -R | \
head -"$PERTYPE" -> ag.nonsplicesite

python ../scripts/nonsplicesites.py ../hg19/hg19.fa GT | sort -k1,1 -k2,2n | \
bedtools intersect -sorted -a inverse.allexons -b stdin | sort -R | \
head -"$PERTYPE"-> gt.nonsplicesite
