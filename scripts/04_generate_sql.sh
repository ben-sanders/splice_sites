#!/bin/bash
#################################
# 04_generate_sql.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Converts the combined data into per-variant SQL statements so they can be 
# added to the database
#################################

# file variables
REFSEQFILE="../refseq_sites/combined.allexons"
RANDFILE="../random_sites/combined.nonsplicesites"
LITFILE="../literature_sites/raw/literature_variants_with_header.txt"

# check the files exist, quit if one doesn't
if [ ! -e "$REFSEQFILE" ]
then
    echo "Cannot open "$REFSEQFILE""
fi

if [ ! -e "$RANDFILE" ]
then
    echo "Cannot open "$RANDFILE""
fi

if [ ! -e "$LITFILE" ]
then
    echo "Cannot open "$LITFILE""
fi

############# Process RefSeq input file #############

# file order for string concat in awk
# 1   2           3       4   5       6       7           8       9       10  11  12
# CHR START_POS   END_POS ID  ACC/DON STRAND  TRANSCRIPT  GENE    MINUS50 REF ALT PLUS50


cat "$REFSEQFILE" | sed s/'\t0\t'/'\tdon\t'/ | sed s/'\t1\t'/'\tacc\t'/ | \
awk 'BEGIN {FS="\\t"; print "BEGIN TRANSACTION;"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, source, gene, type) VALUES ( \47%s\47, \47%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47WILDTYPE\47, \47REFSEQ\47, \47%s\47, \47%s\47);\n", $4, $1, $2, $3, $6, $7, $10, $11, $9, $12, $9$10$12, $9$11$12, 1, $8, $5} END{print"COMMIT;"}' > refseq_variant_upload.sql

rm "$REFSEQFILE"

############# Process random sites #############

# output filename should be random_site_upload.sql
# random sites are formatted the same as RefSeq sites.
cat "$RANDFILE" | sed s/'\t0\t'/'\tdon\t'/ | sed s/'\t1\t'/'\tacc\t'/ | \
awk 'BEGIN {FS="\\t"; print "BEGIN TRANSACTION;"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, source, gene, type) VALUES ( \47%s\47, \47%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47NONSPLICE\47, \47RANDOM\47, \47%s\47, \47%s\47);\n", $4, $1, $2, $3, $6, $7, $10, $11, $9, $12, $9$10$12, $9$11$12, 1, $8, $5} END{print"COMMIT;"}' > random_site_upload.sql

rm "$RANDFILE"

############# Process literature variants #############

# 1     2       3   4       5          	6               7       8       9          	10          11      12          13          14
# ID	Gene	NM_	HGVS	c. number	Affects +-1/2	chr.	strand	position	start_pos	end_pos	upstream	wildtype	mutant

# 15            16      17      18      19
# downstream	Effect	Paper	Method	Comments

grep -v ^ID "$LITFILE" | \
awk 'BEGIN {FS="\\t"; print "BEGIN TRANSACTION;"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, hgvs, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, paper, method, source, gene) VALUES ( \47%s\47, \47chr%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47%s\47, \47%s\47, \47%s\47, \47LITERATURE\47, \47%s\47);\n", $1, $7, $10, $11, $8, $3, $4, $13, $14, $12, $15, $12$13$15, $12$14$15, $6, $16, $17, $18, $2} END{print"COMMIT;"}' > literature_variant_upload.sql

# don't delete the $LITFILE, as it is small, and has to be manually extracted from Excel.

############# Create and update database #############

# delete database if it already exists
rm variantdb.db

# recreate using files created above
touch variantdb.db
sqlite3 variantdb.db '.read create_variantdb.sql'
sqlite3 variantdb.db '.read literature_variant_upload.sql'
sqlite3 variantdb.db '.read refseq_variant_upload.sql'
sqlite3 variantdb.db '.read random_site_upload.sql'
