#
#
#

# this should follow on from 02_get_refseq_sites in the full automated reanalysis pipeline.

# uses the combined.allexons file produced by that script to choose some random non-splice sites as negative controls.

# use bedtools slop to expand each splice site by +-20, to avoid accidentally selecting splice positions.
# NOTE: this will include exons, and so may include ESEs. Should use a list of exons instead? So only intronic sequence is considered? 
# reverse the ranges in combined.allexons, to give a list of all non-splice-site coordinates in the genome

echo "ERROR: This doesn't do anything yet."
