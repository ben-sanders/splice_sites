-- create the empty table for storing the variant database

BEGIN TRANSACTION;
DROP TABLE IF EXISTS variants;

-- db_id is automatically generated  
-- var_id should be how the variant is identified in the original source (i.e. the excel spreadsheet)
CREATE TABLE variants (db_id INTEGER PRIMARY KEY, gene TEXT,
                       var_id TEXT NOT NULL
                       chrom TEXT NOT NULL,
                       pos_start INTEGER NOT NULL,
                       pos_end INTEGER NOT NULL,
                       strand TEXT CHECK(strand IN ('+', '-')),
                       transcript TEXT,
                       hgvs TEXT,
                       ref_allele TEXT,
                       alt_allele TEXT,
                       up_50bp TEXT,
                       down_50bp TEXT,
                       ref_seq TEXT,
                       alt_seq TEXT,
                       canonical_site TEXT,
                       effect TEXT,
                       paper TEXT,
                       method TEXT,
                       comments TEXT,
                       source TEXT);
COMMIT;

/* 
Example queries, for syntax guidance

insert new variant (mandatory minimum shown):

NOTE: CHECK THE REFERENCE USED TO SEE IF CHROMOSOMES SHOULD START WITH CHR
	INSERT INTO variants (var_id, chrom, pos_start, pos_end) VALUES ('001', 'chr1', 12345, 12346)

update variant with new info (reference and alt alleles shown):
    UPDATE variants SET ref_allele='A', alt_allele='G' WHERE var_id='001'
*/

/*
Mappings for variantdb.xls to sqlite update

1   2       3   4       5           6               7       8       9           10          11      12          13          14
ID	Gene	NM_	HGVS	c. number	Affects +-1/2	chr.	strand	position	start_pos	end_pos	upstream	wildtype	mutant

15          16      17      18      19
downstream	Effect	Paper	Method	Comments

don't include comments, as blank fields end up skipping the end of the SQL command.

awk 'BEGIN {FS="\\t"}{printf "INSERT INTO variants (var_id, chrom, pos_start, pos_end, strand, transcript, hgvs, ref_allele, alt_allele, up_50bp, down_50bp, ref_seq, alt_seq, canonical_site, effect, paper, method, comments, source) VALUES ( \47%s\47, \47%s\47, %d, %d, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, \47%s\47, %d, \47%s\47, \47%s\47, \47%s\47, \47LITERATURE\47);\n", $1, $7, $10, $11, $8, $3, $4, $13, $14, 12+13+15, 12+14+15, $6, $16, $17, $18}'

*/

