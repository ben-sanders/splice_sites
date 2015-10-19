###
# convert_to_alamut_format.sh
# author: Ben Sanders
# date:   2015-10-19
###

# Script to convert the excel variant database (exported as tab-separated) into a
# format compatible with input to Alamut batch.
# writes output to console, so can redirect to file.
# does not overwrite input (DON'T DO THAT!) or any other file by default.

# USAGE: ./convert_to_alamut_format.sh [variant_db.txt]

# first check the correct number of inputs
if [ $# != 1 ]
then
    echo "ERROR: expects one input, multiple inputs not supported."
    echo "USAGE: ./convert_to_alamut_format.sh [variant_db.txt]"
    exit
fi

# then check the file exists
if [ ! -f $1 ]
then
    echo "ERROR: specified file does not exist or cannot be opened."
    echo "USAGE: ./convert_to_alamut_format.sh [variant_db.txt]"
    exit
fi

# now we can go through line-by-line and rearrange everything

# NOTE: for testing, starting with one line only.

# input format
# 1	2	3	4	5		6		7	8	9		10		11
# ID	Gene	NM_	HGVS	c. number	Affects +-1/2	chr.	strand	position	upstream	LENGTH CHECK
# 12		13	14		15		16	17	18	19	
# wildtype	mutant	downstream	LENGTH CHECK	Effect	Paper	Method	Comments

# output format
# if using strand, run alamut with --strand 0
# no data shold be blank, not ".", "na", etc...
# user defined fields (cols 8+) are ignored and just reported in output.
# insertions + deletions: - for the version with less sequence (e.g. ALT for deletions, REF for insertions)

#########
# WARNING: It looks like Alamust batch might not support indels?
#########

# 1	2	3	4	5	6		7		8
# ID	CHR	POS	REF	ALT	STRAND(opt)	TRANSCRIPT(opt) USER_DEFINED(opt)

# input - output mapping
# 1  - 1
# 7  - 2
# 9  - 3
# 12 - 4
# 13 - 5
# 8  - 6
# 3  - 7

# This rejigs the format correctly, but there may need to be some further cleaning after?
# e.g. how does is handle multi-base changes? DOES it handle multi-base changes?
# if not, then may need to get rid of positions with dashes (i.e. range-based positions)
cat $1 | awk 'BEGIN{FS="\\t"}{print $1, "\t", $7, "\t", $9, "\t", $12, "\t", $13, "\t", $8, "\t", $3}'
