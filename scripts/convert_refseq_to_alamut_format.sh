###
# convert_refseq_to_alamut_format.sh
# author: Ben Sanders
# date:   2015-10-19
###

# run on allexons.bed
# this file contains the chromosome, position, strand, and identifier for all refseq exon-inton boundaries

# input format
# 1     2           3       4   5       6
# CHR   POS_START   POS_END ID  TYPE    STRAND

# need to get the ref base from somewhere!
# could run getfasta on this, should return >ID\nREF
# then would need to match ID to ID
# better to put into a database?
# could then add by iterating and using ID as primary key.

# output format
# if using strand, run alamut with --strand 0
# no data shold be blank, not ".", "na", etc...
# user defined fields (cols 8+) are ignored and just reported in output.
# insertions + deletions: - for the version with less sequence (e.g. ALT for deletions, REF for insertions)

#########
# WARNING: It looks like Alamust batch might not support indels?
#########

# 1     2	3	4	5	6		    7		        8
# ID    CHR	POS	REF	ALT	STRAND(opt)	TRANSCRIPT(opt) USER_DEFINED(opt)

# input - output mapping
