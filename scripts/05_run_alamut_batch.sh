#!/bin/bash
#################################
# 05_run_alamut_batch.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Extract variant information from the database and run Alamut batch
#################################



# Alamut batch wants data in tab delimited files, but not in BED format.
# From the Alamut documentation, the fields are:

# 1. Variant id (anything)
# 2. Chromosome (1-22, X, Y)
# 3. Genomic position
# 4. Reference nucleotide(s) (ACGT, or '-' for insertions)
# 5. Mutated nucleotide(s) (ACGT, or '-' for deletions)
# 6. Optional strand (1/+ or -1/-), used if --strand parameter is set to 01
# 7. Optional transcript id, used if --spectrans parameter is specified
# 8. Optional user-defined fields (e.g. heterozygosity, number of reads, etc).

# all variant data is in the database/variantdb.db, and sqlite3 database

# to output as tab-separated:
sqlite3 ../database/variantdb.db '.mode tabs'
# then whatever output query is wanted
sqlite3 ../database/variantdb.db '
