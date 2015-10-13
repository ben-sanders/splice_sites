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

1. Download list of RefSeq genes/transcripts from UCSC table browser.
 * group: "Genes and Gene Predictions"
 * track: "RefSeq Genes"
 * table: "refGene"
 * region: "genome"
2. The downloaded file contains a list of RefSeq transcripts, and the genes with which the are associated. The script make_exon_list.sh converst this list into a list of splice-sites extracted from the transcript with the most splice sites from each gene.
 * getrefseqexons.py is run, which finds the longest (most exons) transcript, extracts the exon/intron boundary locations, and writes this information to stdout in BED format.
 * bedtools slop is used to increase the size of the bed region to 50bp either side of the single base given by getrefseqexons.py. This is for flexibility, o the sie of selection can be easily changed. It also has built in controls to ensure a selection cannont exceed chromosome boundaries.
 * bedtools getfasta uses the BED file to extract the genomic sequence from a given reference file (**NOTE:** *hg19 reference file is too big to include in the repo.*). FASTA identifier lines are automatically generated
 * Since the FASTA is extracted from a single reference, any genes on the reverse strand will be incorrectly oriented. strandrevcomp.py reverse complements these sequences so that they are correct for the transcript orientation.
3. Generate random / decoy sequences
 * Uses the Python random module to select random bases to create sequences of a defined length. The number of sequences is also specifiable, but should match the number of refseq splice sites. Ideally, these sequences will not contain sequences matching splice sites, but may also give an indication of the occurence of high-scoring motifs by chance. Real functional sites will be dependent on other features that should be missing. Odds of any one n-base sequence are 4^n, so longer sequences are less likely to occur by chance. However, as there is a degree of flexibility in splice-sites, the actual odds of a site occuring by chance are lower than this. (If there was no flexibility, this would be a much easier task!)
4. Run predictors
 * Before doing this, check that no new predictors have been published
 * Either:
  * Use a custom script to run each tool, so that the output format is controllable, or:
  * Use alamut batch for the majority, and make a script to parse and interpret this output
5. Interpret results
 * Identifying significant changes:
  * Best practice guidelines suggest >10% change, presumable 10% reduction!
  * probably needs to be considered along with other sites - 10% change might still be the highest score!
   * <10% reduction == No effect
   * >10% reduction == Effect
   * Use this arbitrary cutoff or attempt to learn best params?
  * Consider position of true splice site in list of scores (for tools that scan and score, e.g. mespy)
   * wt highest > var highest == No effect
   * wt highest > var NOT highest == Effect
   * wt NOT highest > var NOT highest == Possible effect depending on order + magnitude of score change.
   * Use relative score (i.e. difference between wt and var rank) as a modifier somehow?
    * e.g. >10% predicitve score reduction AND < -1 position score == Effect
    * score reduction but no position change == probably no effect
  * For literature variants, sequence is centered around variant, not splice site. This may cause some problems. Might be good to establish the position of the actual splice site, so we can look at its score + ranking.
  * Is it ok to assume that the highest scoring position in the wildtype sequence is the splice site? Probably not, but how best to account for the position?

