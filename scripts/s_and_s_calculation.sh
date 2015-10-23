
# Uses the RefGene list to get all splice sites, then splits by type
# Selects the appropriate ranges for each type to match the Shapiro and Senepath (1987) splice-site matrices
# then runs a Python script to work out the distribution of bases at each position.

# USAGE: ./s_and_S_calculation <refgenelist>

python ~/Documents/Ben_project/splice_sites/scripts/getrefseqexons.py "$1" > refseq_splicesites.bed

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

# split the records by the type of site
# need to correct for orientation - e.g. reverse complement all -ve strand sites, so separate them too.

cat plusminus1.allexons.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 0 && $6 == "+") print}' > positive.donors.bed
cat plusminus1.allexons.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 0 && $6 == "-") print $1 "\t" $2 "\t" $3 "\t" $4 "-\t" $5 "\t" $6 "\t" $7 "\t" $8}' > negative.donors.bed
cat plusminus1.allexons.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 1 && $6 == "+") print}' > positive.acceptors.bed
cat plusminus1.allexons.bed | awk 'BEGIN{FS="\\t"}{if ($5 == 1 && $6 == "-") print $1 "\t" $2 "\t" $3 "\t" $4 "-\t" $5 "\t" $6 "\t" $7 "\t" $8}' > negative.acceptors.bed

# Donors
# 3 bases of exon, 6 bases of exon

bedtools slop -i positive.donors.bed -g /home/genetics/Documents/ReferenceGenome/hg19.fa.fai -l 4 -r 5 > positive.ss.donors.bed
bedtools slop -i negative.donors.bed -g /home/genetics/Documents/ReferenceGenome/hg19.fa.fai -l 6 -r 3 > negative.ss.donors.bed

# Acceptors
# 13 bases of intron, 1 of exon
bedtools slop -i positive.acceptors.bed -g /home/genetics/Documents/ReferenceGenome/hg19.fa.fai -l 14 -r 1 > positive.ss.acceptors.bed
bedtools slop -i negative.acceptors.bed -g /home/genetics/Documents/ReferenceGenome/hg19.fa.fai -l 2 -r 13 > negative.ss.acceptors.bed
    
# retreive FASTA sequences
bedtools getfasta -fi /home/genetics/Documents/ReferenceGenome/hg19.fa -name -bed positive.ss.donors.bed -fo positive.ss.donors.fa
bedtools getfasta -fi /home/genetics/Documents/ReferenceGenome/hg19.fa -name -bed negative.ss.donors.bed -fo negative.ss.donors.fa
bedtools getfasta -fi /home/genetics/Documents/ReferenceGenome/hg19.fa -name -bed positive.ss.acceptors.bed -fo positive.ss.acceptors.fa
bedtools getfasta -fi /home/genetics/Documents/ReferenceGenome/hg19.fa -name -bed negative.ss.acceptors.bed -fo negative.ss.acceptors.fa

# reverse complement the negatives
# this script prints some stuff to stderr, so get rid of it
python ../scripts/strandrevcomp.py negative.ss.donors.fa > revcomp.negative.ss.donors.fa 2> /dev/null
python ../scripts/strandrevcomp.py negative.ss.acceptors.fa > revcomp.negative.ss.acceptors.fa 2> /dev/null

# combine into single files, without and IDs

grep -v "^>" positive.ss.donors.fa > combined.donors.intermediate
grep -v "^>" revcomp.negative.ss.donors.fa >> combined.donors.intermediate
grep -v "^>" positive.ss.acceptors.fa > combined.acceptors.intermediate
grep -v "^>" revcomp.negative.ss.donors.fa >> combined.acceptors.intermediate

# there are a few lowercase sequences
tr '[:lower:]' '[:upper:]' < combined.donors.intermediate > combined.donors
tr '[:lower:]' '[:upper:]' < combined.acceptors.intermediate > combined.acceptors

# clean up BED files and unprocess FASTA
rm *.bed *.fa *.intermediate

python ../scripts/s_and_s_pwm.py combined.donors
python ../scripts/s_and_s_pwm.py combined.acceptors


