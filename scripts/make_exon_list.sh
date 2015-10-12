# from a given RefSeq gene list, extract splice sites BED file
# the gene list location if defined in the python script, but it probably
# doesn't need to be
python ~/Documents/misc/splice_sites/scripts/getrefseqexons.py > refseq_splicesites.bed

# bedtools slop
bedtools slop \
    -i refseq_splicesites.bed \
    -g ~/Documents/misc/splice_sites/reference_genome/hg19/hg19.genome \
    -b 50 > slop.refseq_splicesites.bed

# bedtools getfasta
bedtools getfasta \
    -name \
    -fi ~/Documents/misc/splice_sites/reference_genome/hg19/hg19.fa \
    -bed slop.refseq_splicesites.bed \
    -fo fasta.refseq_splicesites.fa

# reverse complement negative strand sequences
python ~/Documents/misc/splice_sites/scripts/strandrevcomp.py > revcomp.refseq_splicesites.fa
