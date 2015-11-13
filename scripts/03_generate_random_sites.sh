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
GENOMEFILE="../hg19/hg19.genome"
GENOMEFASTA="../hg19/hg19.fa"
#NOTE: This takes a long time to run, as it has to iterate over the entire
#      human genome several times

# bedtools complement reverses the selection - returns the regions of the
# genome file that aren't covered by intervals in the combined.allexons bedfile.
cat "$EXONFILE" | bedtools slop -i stdin -g "$GENOMEFILE" -b 20 | \
bedtools complement -i stdin -g "$GENOMEFILE" > inverse.allexons

# Generate bed file of all AG/GT sites in the genome, get rid of
# all that intersect with slopped combined.allexons ranges, and select at random
# ensure a matching split between donors and acceptors (e.g. half of each,
# although this could be worked out exactly if that seems inaccurate - but how 
# could there be introns with donors but not acceptors? Since I removed first
# and last exon boundaries...)

# this is slow, so check if either file already exists first.
if [ ! -e combined.nonsplicesites ]
then
    python ../scripts/nonsplicesites.py "$GENOMEFASTA" AG | \
    bedtools intersect -a stdin -b inverse.allexons | sort -R | \
    head -"$PERTYPE" > nonsplicesites.bed
    python ../scripts/nonsplicesites.py "$GENOMEFASTA" GT | \
    bedtools intersect -a stdin -b inverse.allexons | sort -R | \
    head -"$PERTYPE" >> nonsplicesites.bed
fi

# This outputs: CHR START_POS   END_POS
# start and end are flanking the AG/GT dinucloetide
# to match the refseq sites, want to include only the +- 1 position
# fields wanted by database
# var_id, chrom, pos_start, pos_end, strand, transcript, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, source, gene, type
# to get this copy + modify code from 02_get_refseq_sites

# replace X,Y, and M with numbers so the sort correctly. 
# M should come first in hg19, so set to 0 
sed -i s/"chrM"/"chr0"/ nonsplicesites.bed
sed -i s/"chrX"/"chr23"/ nonsplicesites.bed
sed -i s/"chrY"/"chr24"/ nonsplicesites.bed

# do the sort here, then reset chromosome names
# sort by chromosome and start positions (columns 1 and 2)
sort nonsplicesites.bed  -k1.4n,2n  > sorted.nonsplicesites.bed 

sed -i s/"chr0"/"chrM"/ sorted.nonsplicesites.bed 
sed -i s/"chr23"/"chrX"/ sorted.nonsplicesites.bed 
sed -i s/"chr24"/"chrY"/ sorted.nonsplicesites.bed 

cat sorted.nonsplicesites.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 0 && $6 == "+" ) print $1 "\t" $2+1 "\t" $3+1 "\t" $4 "\t" $5 "\t" $6 "\tTSCRIPT_PLACEHOLDER\tGENE_PLACEHOLDER"; else if ($5 == 1 && $6 == "-") print $1 "\t" $2+1 "\t" $3+1 "\t" $4 "\t" $5 "\t" $6 "\tTSCRIPT_PLACEHOLDER\tGENE_PLACEHOLDER"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tTSCRIPT_PLACEHOLDER\tGENE_PLACEHOLDER"}' > plusminus1.nonsplicesites.bed

#The +-1 target base is then covered by removing 1 from the start index. This could have been included above, but 

cat plusminus1.nonsplicesites.bed | awk 'BEGIN{FS="\\t"}{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8;}' > wt.nonsplicesites.bed

# -name adds the ID from the bed file as the fasta ID.	
#check the ordering is the same by comparing with id column from bed file.

#NOTE: the minus50 bed file, when converted to FASTA, misses the base between the end of the range and the splice-site base. The end position of the minus50 needs 1 adding

bedtools flank -i plusminus1.nonsplicesites.bed -g "$GENOMEFILE" -l 50 -r 0 > minus50.nonsplicesites.bed
bedtools flank -i plusminus1.nonsplicesites.bed -g "$GENOMEFILE" -l 0 -r 49 | awk 'BEGIN{FS="\\t"}{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}'> plus50.nonsplicesites.bed

##### will need to replace this with genome fasta later ######
bedtools getfasta -fi "$GENOMEFASTA" -name -bed wt.nonsplicesites.bed -fo wt.nonsplicesites.fa
bedtools getfasta -fi "$GENOMEFASTA" -name -bed minus50.nonsplicesites.bed -fo minus50.nonsplicesites.fa
bedtools getfasta -fi "$GENOMEFASTA" -name -bed plus50.nonsplicesites.bed -fo plus50.nonsplicesites.fa

# remove the id lines from fasta outputs with grep
#	
grep -v "^>" wt.nonsplicesites.fa > wt.nonsplicesites
grep -v "^>" minus50.nonsplicesites.fa > minus50.nonsplicesites
grep -v "^>" plus50.nonsplicesites.fa > plus50.nonsplicesites

# NOTE: All fasta output is from the +ve strand. Sequences may need to be reversed to get the true sequence for the gene.
# NOTE: Because of the way the reference genome is formatted, some bases are lowercase (? in repeated regions?). Convert to uppercase with:

tr '[:lower:]' '[:upper:]' < wt.nonsplicesites > upper.wt.nonsplicesites
tr '[:lower:]' '[:upper:]' < minus50.nonsplicesites > upper.minus50.nonsplicesites
tr '[:lower:]' '[:upper:]' < plus50.nonsplicesites > upper.plus50.nonsplicesites

# Create an alternative allele for use in alamut. Since A or T should abolish splicing in the majority of cases, use them. If alread A or T, swap for the other.

cat upper.wt.nonsplicesites | awk '{if ($1 == "C" || $1 == "T") print "A"; else print "T"}' > alt.nonsplicesites

# Combine into one output:
paste plusminus1.nonsplicesites.bed upper.minus50.nonsplicesites  upper.wt.nonsplicesites  alt.nonsplicesites  upper.plus50.nonsplicesites  > combined.nonsplicesites

# cleanup
rm *.bed *.fa *wt.nonsplicesites *50.nonsplicesites alt.nonsplicesites
