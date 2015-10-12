"""

Wrapper for the mespy.py script, to enable running on every sequence in a file
comparing 2 sequences at a go - one for wild-type, one for mutant/non-splice.
Will do a comparison between the two results - >10% reduction in strength, or 
another site is not stronger

"""

import mespy
from operator import itemgetter

test = "CCGCATTTGTACTGTGTTAAGCAATAAAAAGTATACCTTTTTGTTGTCAGGCTTGAGAAAGAGATCCAAGATCTTGAAAAAGCTGAACTGCAAATCTCAA"
test2 = "ACTACGAGCATCGACGGCGACGACTATCTACTACTACTAGCAGCTACATTATATAGCAGCACTATCAGCTCTAGCTGCGATGCTAGCTACGATCGTATAT"

scores = mespy.window(test)
scores = sorted(scores, key=itemgetter(3), reverse=True)

scores2 = mespy.window(test2)
scores2 = sorted(scores2, key=itemgetter(3), reverse=True)
    
# check the highest scoring
print "wt: %.3f, rand: %.3f" % (scores[0][3], scores2[0][3])


# check the differences by position
