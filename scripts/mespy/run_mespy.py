"""

Wrapper for the mespy.py script, to enable running on every sequence in a file
comparing 2 sequences at a go - one for wild-type, one for mutant/non-splice.
Will do a comparison between the two results - >10% reduction in strength, or 
another site is not stronger

Score site position by change in ranking, e.g. drop one place = -1, stay the same = 0, improve 2 = 2, etc...

"""

import mespy
from operator import itemgetter

def score_ranking(wt_site, wt_scores, var_scores):
    """ 
    given the wt site (which may just be the highest scoring), locate in the variant list, and
    return the position difference between the two lists.
    """
    print wt_site
    print "---------------------"
    wt_rank = find_site(wt_site, wt_scores)
    var_rank = find_site(wt_site, var_scores)
    print "rank score:", wt_rank-var_rank

    wt_score = wt_site[3]
    var_score = var_scores[var_rank][3]
    print wt_score, var_score
    # take wt from var, so increase == positive, decrease == negative
    score_change = var_score - wt_score
    print score_change
    print var_score / wt_score * 100

def find_site(test_site, score_list):
    """
    returns the position of a site within a given list of sites, ranked by score.
    Need to search by position, not the whole record - won't be found in var list!
    If not found, for whatever reason (there shouldn't be a length restriction to the list), returns -1.
    """
    test_site_position = test_site[1]
    position = 0
    for site in [x[1] for x in score_list]:
        if site == test_site_position:
            return position
        else:
            position += 1
    return -1

def generate_scores(seq1, seq2):
    """
    Runs mespy window across each sequence, and returns lists ranked by score (from highest to lowest)
    """    
    scores = mespy.window(test)
    scores = sorted(scores, key=itemgetter(3), reverse=True)

    scores2 = mespy.window(test2)
    scores2 = sorted(scores2, key=itemgetter(3), reverse=True)
    
    return scores, scores2

test = "CCGCATTTGTACTGTGTTAAGCAATAAAAAGTATACCTTTTTGTTGTCAGGCTTGAGAAAGAGATCCAAGATCTTGAAAAAGCTGAACTGCAAATCTCAA"
test2 = "ACTACGAGCATCGACGGCGACGACTATCTACTACTACTAGCAGCTACATTATATAGCAGCACTATCAGCTCTAGCTGCGATGCTAGCTACGATCGTATAT"

wt_scores, var_scores = generate_scores(test, test2)
score_ranking(wt_scores[0], wt_scores, var_scores)


# check the differences by position
