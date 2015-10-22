###
# run_alamut.sh
# author: Ben Sanders
# date:   2015-10-19
###

# Use this file to run Alamut, in order to ensure the same settings are used every time.
# Could connect with convert_db_to_alamut_format.sh, but once converted don't need to run that every time!

# USAGE: ./run_alamut [varfile] [runID]
# NOTE: runID will be checked and should be unique

# outputs: annotation file, failed variants file, log file.

# append runID to file names for output files
varfile="$1"
annfile="$2_annotation"
failfile="$2_failed_variants"

# assemble the alamut command here
run_alamut="alamut-batch --in $varfile --ann $annfile --unann $failfile --strand 0 --spectrans --ignoreInputErrors --nomispred"
# --in                  = input file
# --ann                 = annotated output file
# --unann               = failed variants - could not be annotated
# --strand 0            = strand set per variant
# --spectrans           = use the specific transcript given for each variant
# --ignoreInputErrors   = runs if some inputs are incorrect - ?writes to failed file?
# --nomispred           = disable missense predictions (only interesting in splicing)

# first check the correct number of inputs
if [ $# != 2 ]
then
    echo "ERROR: expects two inputs."
    echo "USAGE: ./run_alamut [varfile] [runID]"
    exit
fi

# check the input file can be opened
if [ ! -f $1 ]
then
    echo "ERROR: specified input variant file does not exist or cannot be opened."
    echo "USAGE: ./run_alamut [varfile] [runID]"
    exit
fi

# then check the output files don't exist - don't want to overwrite anything
if [ -f $annfile ] || [ -f $failfile ]
then	
    echo "ERROR: Output files already exist, must not overwrite existing files. Please repeat with new runID"
    echo "USAGE: ./run_alamut [varfile] [runID]"
    exit
fi

echo "all preconditions passed, this is where I'd run Alamut..."
echo $run_alamut
