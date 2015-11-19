# try getting coordinates and taking even more flanking sequence with bedtools
# use -s in getfasta to get reverse complement for -ve strand sequences
sqlite3 ../database/variantdb.db 'SELECT chrom, pos_start, pos_end, var_id, db_id, strand FROM variants WHERE source="REFSEQ";' | \
sed s/"|"/"\t"/g | \
bedtools slop -i stdin -g ../hg19/hg19.genome -b 500 | \
bedtools getfasta -fi ../hg19/hg19.fa -bed stdin -fo stdout -name -s > genesplicer.refseq

./../scripts/genesplicer.sh genesplicer.refseq > out.genesplicer.refseq


cd ../random_sites
# for now, select a random subset
sqlite3 ../database/variantdb.db 'select chrom, pos_start, pos_end, var_id, db_id, strand from variants where source="RANDOM";' | \
sed s/"|"/"\t"/g | \
bedtools slop -i stdin -g ../hg19/hg19.genome -b 500 | \
bedtools getfasta -fi ../hg19/hg19.fa -bed stdin -fo stdout -name -s > genesplicer.random

./../scripts/genesplicer.sh genesplicer.random > out.genesplicer.random
