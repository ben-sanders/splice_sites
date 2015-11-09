# Convert the combined data into per-variant SQL statements so they can be added to the database

cat combined.allexons | awk 'BEGIN {FS="\\t"; print "BEGIN TRANSACTION;"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, source, gene) VALUES ( \47%s\47, \47%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47WILDTYPE\47, \47REFSEQ\47, \47%s\47);\n", 
$4, $1, $2, $3, $6, $7, $10, $11, $9, $12, $9$10$12, $9$11$12, 1, $8} END{print"COMMIT;"}' > refseq_variant_upload.sql

# Also do the literature variants

grep -v ^ID literature_variants_with_header.txt | awk 'BEGIN {FS="\\t"; print "BEGIN TRANSACTION;"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, hgvs, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, paper, method, source, gene) VALUES ( \47%s\47, \47chr%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47%s\47, \47%s\47, \47%s\47, \47LITERATURE\47, \47%s\47);\n", $1, $7, $10, $11, $8, $3, $4, $13, $14, $12, $15, $12$13$15, $12$14$15, $6, $16, $17, $18, $2} END{print"COMMIT;"}' > literature_variant_upload.sql

# then create + update database

# Not sure if this can be done automatically?

# $ sqlte3 variantdb.db
# sqlite3> .read create_variant_db.sql
# sqlite3> .read literature_variant_upload.sql
# sqlite3> .read refseq_variant_upload.sql
