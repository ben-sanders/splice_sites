"""

Wrapper for the mespy.py script, to enable running on every sequence in a file
comparing 2 sequences at a go - one for wild-type, one for mutant/non-splice.
Will do a comparison between the two results - >10% reduction in strength, or 
another site is not stronger

Score site position by change in ranking, e.g. drop one place = -1, stay the same = 0, improve 2 = 2, etc...

"""

import mespy
from operator import itemgetter
import sys

def score_site(wt_site, wt_scores, var_scores):
    """ 
    given the wt site (which may just be the highest scoring), locate in the variant list, and
    return the position difference between the two lists.

    Inputs:
        wt_site    - site (ideally true wild-type splice site, being tested)
        wt_scores  - list of scores from wild-type sequence
        var_scores - list of scores from variant sequence 
    Outputs:    
        (rank change, percentage change)
    """
    
    # find the rank (position) of the test site in wt and variant score lists
    wt_rank = find_site(wt_site, wt_scores)
    var_rank = find_site(wt_site, var_scores)
    
    # use rank to find the exact score (rank is 0 indexed, as are the score lists)
    wt_score = wt_site[3]
    var_score = var_scores[var_rank][3]
    # take wt from var, so increase == positive, decrease == negative
    score_diff = var_score - wt_score
    score_change = (score_diff / abs(wt_score)) * 100
    rank_change = wt_rank - var_rank
    return rank_change, score_change, wt_score, var_score

def find_site(test_site, score_list):
    """
    returns the position of a site within a given list of sites, ranked by score.
    Need to search by position, not the whole record - won't be found in var list!
    If not found, for whatever reason (there shouldn't be a length restriction to the list), returns -1.
    """
    # get position of the test splice site, this will allow comparison to the same point
    test_site_position = test_site[1]
    position = 0
    
    # count through the list until you find the position you're looking for
    # then return the index it is at.
    for site in [x[1] for x in score_list]:
        if site == test_site_position:
            return position
        else:
            position += 1
    # if not found, return -1
    return -1

def generate_scores(seq1, seq2, mode):
    """
    Runs mespy window across each sequence, and returns lists ranked by score (from highest to lowest)
    window output is (TYPE, POS_IN_SEQ, SEQ_AS_EVALUATED, SCORE)
    """    
    scores = mespy.window(seq1, mode)
    scores = sorted(scores, key=itemgetter(3), reverse=True)

    scores2 = mespy.window(seq2, mode)
    scores2 = sorted(scores2, key=itemgetter(3), reverse=True)
    
    return scores, scores2

def run_mespy(seq1, seq2, mode="both"):
    """
    Just a wrapper for testing purposes.
    """
    wt_scores, var_scores = generate_scores(seq1, seq2, mode)
    numwt = len(wt_scores)
    numvar = len(var_scores)
    if numwt != 0:
        rank_change, score_change, wt_score, var_score = score_site(wt_scores[0], wt_scores, var_scores)
        position = wt_scores[0][1]
        wt_score = wt_scores[0][3]
        return position, rank_change, score_change, wt_score, var_score
    else:
        return -1, -999, -999
        
def read_input(fname):
    """
    Run mespy on variants in file from variant database
    Doesn't support vcf, as it needs to include the full sequence.
    Input is in this format:
    
    DB_ID   SITE_ID REF_SEQ ALT_SEQ TYPE
    
    output (to stdout) is:
    
    DB_ID   SITE_ID POS RANK_CHANGE  SCORE_CHANGE(%)
    """
    print "DB_ID\tSITE_ID\tPOS\tRANK_CHANGE\tWT\tVAR\tSCORE_CHANGE"
    f = open(fname, "r")
    for line in f:
        line = line.strip().split()
        ref = line[2]
        alt = line[3]
        if line[4] == "acc":
            # look at 3' acceptor motif only
            mode = 3
        elif line[4] == "don":
            # look at 5' donor motif only
            mode = 5
        else:
            # site type not specified, so examine both
            mode = "both"
        
        if "N" in line[2] or "N" in line[3]:
            # N characters cause the scoring to crash
            print "%s\t%s\tN\tN\tN\tN\tN"
        else:        
            position, rank_change, score_change, wt_score, var_score = run_mespy(ref, alt, mode)
            print "%s\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f" % (line[0], line[1], position, rank_change, wt_score, var_score, score_change)    
    
if __name__ == "__main__":
    #test1 = "CCCGCGGACCCTCCCTCCCGGCCTTCCGCCACCGGCGCGGGCGCAACTCACCGGGCATCAGCTCTTCCGGCTCCCTCATGCCACGGGCAGTACGGGCAGCC"
    #test2 = "CCCGCGGACCCTCCCTCCCGGCCTTCCGCCACCGGCGCGGGCGGGGCTCAACGGGCATCAGCTCTTCCGGCTCCCTCATGCCACGGGCAGTACGGGCAGCC"
    #mode = 3
    print sys.path[0]
    read_input(sys.argv[1])
    #run_mespy(test1, test2, mode)
