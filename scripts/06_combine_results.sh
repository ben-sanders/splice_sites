#!/bin/bash
#################################
# 05_combine_results.sh
# ------------------
# author: Ben Sanders
# ben.sanders@salisbury.nhs.uk
# ------------------
# 2015-11-10
# ------------------
# Integrates the output from the various analysis tools (Alamut, Spidex, etc),
# works out the changes between wt and variant, and runs scripts to compare +
# analyse the whole dataset.
#################################

#TODO

# will probably hacve to spin out analysis scripts for every individual tool,
# prior to combining results in here 
# e.g. run alamut, spidex, etc. then join results on id
# If use the db ID wherever there's an identifier:
#   1. it won't get truncated - SPIDEX limits IDs to 20 characters
#   2. it's guaranteed unique for each variant (created as added to db)
#   3. easily sortable (+ checkable)
#   4. will sort exactly the same, so can join columns with cut and paste


# alamut output is quite long, with a large number of fields
# presumably it keeps them if there is no data? i.e. a blank column rather than
# skipping the entire column?

#  1. id
#  2. chrom
#  3. position (from input)
#  4. reason for failure
#  5. gene
#  6. gene ID (HGNC)
#  7. transcript
#  8. strand
#  9. transcript length
# 10. protein ID
# 11. uniprot id
# 12. variant type
# 13. coding effect (synonymous, missense, etc.)
# 14. variant location (relative to gene - upstream, utr, splice site, etc.)
# 15. genome assembly
# 16. gDNA start
# 17. gDNA end
# 18. HGVS genome nomenclature (for variant)
# 19. cDNA start
# 20. cDNA end
# 21. HGVS cDNA nomenclature
# 22. HGVS protein nomenclature
# 23. alternative protein nomenclature (full for synonymous, rather than p.=)
# 24. exon (numbered? nearest if intronic variant)
# 25. intron (see above!)
# 26. OMIM id (for gene, if present?)
# 27. distance to nearest splice site
# 28. type of nearest splice site
#############################################################
# These are the fields of interest - splice predictor results
#############################################################
# 29. WT SSFL score
# 30. WT MES score
# 31. WT NNSPLICE score
# 32. WT GeneSplicer score
# 33. WT HSF score
# 34. var SSFL score
# 35. var MES score
# 36. var NNSPLICE score
# 37. var GeneSplicer score
# 38. var HSF score
# 39. predicted change to nearest splice site
# 40. predicted effect at variant site
##############################################################
# 41. protein domain 1 ???
# 42. protein domain 2 ???
# ...
# ...
# There are, frankly, a ridiculous number of output fields. I can't be bothered
# to transcribe them all here, but the useful ones are highlighted.
################################################################################ 
