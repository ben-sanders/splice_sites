# Splice Site Predictor Evaluation Project

## directory structure

### literature_sites
Contains the list of splice-site variants derived from literature sources. All have been tested at RNA level for an accurate characteristation of their impact on splicing.

**NOTE:** *Doesn't actually have anything in it yet...*

### random_sites
For use with refseq_sites, contains a matched length list of randomly generated DNA sequences that, ideally, should not contain functional splice motifs.

**NOTE:** *Doesn't actually have anything in it yet...*

### refseq_sites
Contains a list of all splice sites across the human genome, derived from RefSeq records 

#### exons
contains the BED files (at various stages of processsing) needed to extract the splice-site sequences.

#### gene_list
Original list of transcripts (more than one for most genes, later processing selected one transcript per gene with the most exons) downloaded from USCS table browser.
Exact search parameters in README file?

### scripts
Storage for the scripts used in the analysis processs.

### test
Temporary folder for now, used while testing and fine-tuning the site extraction script!

## Analysis process
Fully detailed in the writeup... Which should probably have its own folder too!

1. ...
2. ...
