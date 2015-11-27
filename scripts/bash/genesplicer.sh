# 5. GeneSplicer
# --------------

# GeneSplicer takes fast inputs, but doesn't seem to account for multiple
# samples - it concatenates them all together!?
# need to feed it individual files, which obviously would be a pain so will
# extract database as one file and then pull out the individual sequences for 
# GeneSplicer

# will work this out in a separate script, then copy in here.

# extract from database, create two fasta records - wt and var, using the db
# and variant ids. Add the strand to the end of the id, so that strandrevcomp.py
# can be used to reverse complement negative sequences (to get the sequence 
# relative to the actual gene)

# NOTE: strandrevcompy.py doesn't work on streams, so have to write to an
#       intermediate file.

# This is where GeneSplicer will actually run
echo -e "ID\tSCORE"
while read p; do
    if [ "${p:0:1}" == ">" ]
    then
        varid="$p"
    else
        varfasta="$p"
        echo -e "$varid""\n""$varfasta" > temp.genesplicer
        # run genesplicer, and use awk to select the prediction nearest the 
        # centre point (i.e. where we know the true splice-site is)
        gsscore=`genesplicer temp.genesplicer ../../project_tools/GeneSplicer/human/ -a -1000 -d -1000 2> /dev/null | \
        awk '{if ($1 > 490 && $2 <510) print $3}'`
        # GeneScanner doesn't detect all sites, so check if there is a score
        # before printing. As for MES, give an unpredicted site 999, well out of
        # the output range, and can be run through the same R script.
        if [ "$gsscore" == "" ]
        then
            echo -e ${varid:1}"\t999"
        else
            echo -e ${varid:1}"\t""${gsscore%%$'\n'*}"
        fi
    fi   
done < "$1" 

