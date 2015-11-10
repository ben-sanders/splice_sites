#!/bin/bash
#################################
# 02_get_refseq_sites.sh
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

python ../scripts/getrefseqexons.py "$1" > refseq_splicesites.bed

# replace X,Y, and M with numbers so the sort correctly. 
# M should come first in hg19, so set to 0 
sed -i s/"chrM"/"chr0"/ refseq_splicesites.bed 
sed -i s/"chrX"/"chr23"/ refseq_splicesites.bed 
sed -i s/"chrY"/"chr24"/ refseq_splicesites.bed 

# do the sort here, then reset chromosome names
# sort by chromosome and start positions (columns 1 and 2)
sort refseq_splicesites.bed -k1.4n,2n  > sorted.refseq_splicesites.bed

sed -i s/"chr0"/"chrM"/ sorted.refseq_splicesites.bed 
sed -i s/"chr23"/"chrX"/ sorted.refseq_splicesites.bed 
sed -i s/"chr24"/"chrY"/ sorted.refseq_splicesites.bed 

cat sorted.refseq_splicesites.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 0 && $6 == "+" ) print $1 "\t" $2+1 "\t" $3+1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8; else if ($5 == 1 && $6 == "-") print $1 "\t" $2+1 "\t" $3+1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8; else print}' > plusminus1.allexons.bed

#The +-1 target base is then covered by removing 1 from the start index. This could have been included above, but 

cat plusminus1.allexons.bed | awk 'BEGIN{FS="\\t"}{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8;}' > wt.allexons.bed

# -name adds the ID from the bed file as the fasta ID.	
#check the ordering is the same by comparing with id column from bed file.

#NOTE: the minus50 bed file, when converted to FASTA, misses the base between the end of the range and the splice-site base. The end position of the minus50 needs 1 adding

bedtools flank -i plusminus1.allexons.bed -g ../hg19/hg19.genome -l 50 -r 0 > minus50.allexons.bed
bedtools flank -i plusminus1.allexons.bed -g ../hg19/hg19.genome -l 0 -r 49 | awk 'BEGIN{FS="\\t"}{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}'> plus50.allexons.bed

bedtools getfasta -fi ../hg19/hg19.fa -name -bed wt.allexons.bed -fo wt.allexons.fa
bedtools getfasta -fi ../hg19/hg19.fa -name -bed minus50.allexons.bed -fo minus50.allexons.fa
bedtools getfasta -fi ../hg19/hg19.fa -name -bed plus50.allexons.bed -fo plus50.allexons.fa

# remove the id lines from fasta outputs with grep
	
grep -v "^>" wt.allexons.fa > wt.allexons
grep -v "^>" minus50.allexons.fa > minus50.allexons
grep -v "^>" plus50.allexons.fa > plus50.allexons

# NOTE: All fasta output is from the +ve strand. Sequences may need to be reversed to get the true sequence for the gene.
# NOTE: Because of the way the reference genome is formatted, some bases are lowercase (? in repeated regions?). Convert to uppercase with:

tr '[:lower:]' '[:upper:]' < wt.allexons > upper.wt.allexons
tr '[:lower:]' '[:upper:]' < minus50.allexons > upper.minus50.allexons
tr '[:lower:]' '[:upper:]' < plus50.allexons > upper.plus50.allexons

# Create an alternative allele for use in alamut. Since A or T should abolish splicing in the majority of cases, use them. If alread A or T, swap for the other.

cat upper.wt.allexons | awk '{if ($1 == "C" || $1 == "T") print "A"; else print "T"}' > alt.allexons

# Combine into one output:
paste plusminus1.allexons.bed upper.minus50.allexons upper.wt.allexons alt.allexons upper.plus50.allexons > combined.allexons

# remove intermediary files
rm *.bed *.fa *wt.allexons *50.allexons alt.allexons
