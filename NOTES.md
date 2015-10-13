# For storing notes and comments that don't really fit anywhere else.

* What is the format for Alamut batch input?
 * tab delimited
  * Var ID, chrom, pos, ref, alt, strand (opt.), transcript(opt.), user defined fields (opt.)
 * what does it do if alt is the same? i.e. no variant, I just want the wt scores?
 * may be necessary to just alter each base anyway, and ignore any changes?
 * Allows preferred transcript to be specified - see manual
  * definitely a good idea for the wt sites, make sure it's using the same transcript.

* Also need to know the output format
 * and any additional information that might contain!
 * see manual!
 * important fields:
  * 1) ID    -|
  * 2) chrom  |- for checking against input
  * 3) pos   -|
  * 16) gDNA start  -|- also for checking
  * 17) gDNA end    -|
  * 27) distance to nearest splice site (includes orientation?)
  * 28) wt SSFL score-----------------------|
  * 29) wt MES score (check against mespy)  |
  * 30) wt NNSPLICE score                   |- wild-type scores for nearest splice site
  * 31) wt GeneSplicer score                |
  * 32) wt HSF score -----------------------|
  * 33) var SSFL score--------|
  * 34) var MES score         | 
  * 35) var NNSPLICE score    |- variant sequence scores for nearest splice site
  * 36) var GeneSplice score  |
  * 37) var HSF score --------|
  * 38) predicted change to nearest splice site
* Then, make sure the other tools are going to outouput in a compatible format - e.g. same number of lines per var.
