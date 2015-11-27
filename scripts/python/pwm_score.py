def reversecomplement(sequence):
    """
    Creates the reverse complement of a DNA sequence
    Supports the standard IUPAC nucleotide code (excludes modified bases).
    If an unknown base is seen, rev. comp. sequence will include 'x'. 
    """
    
    # IUPAC complimentary nucleotides. Dictionary used for speed, includes
    # entries for UPPER and lower case.
    IUPAC_NUC = {"A": "T", "a": "t", # A -> T
                 "C": "G", "c": "g", # C -> G
                 "G": "C", "g": "c", # G -> T
                 "T": "A", "t": "a", # T -> A
                 "R": "Y", "r": "y", # A or G -> T or C
                 "Y": "R", "y": "r", # C or T -> G or A
                 "S": "S", "s": "s", # G or C -> G or C
                 "W": "W", "w": "w", # A or T -> A or T
                 "K": "M", "k": "m", # G or T -> A or C
                 "M": "K", "m": "k", # A or C -> G or T
                 "B": "V", "b": "v", # C or G or T -> G or C or A
                 "V": "B", "v": "b", # G or C or A -> C or G or T
                 "D": "H", "d": "h", # A or G or T -> T or C or A
                 "H": "D", "h": "d", # T or C or A -> A or G or T
                 "N": "N", "n": "n", # any base
                 "-": "-"}           # gap
    revcomp = []
    # compliment the sequence
    for base in sequence:
        # get the complimentary code, if one does not exist add 'x'
        revcomp.append(IUPAC_NUC.get(base, "x"))
    # reverse it
    revcomp.reverse()
    # return as a string rather than a list
    return ''.join(revcomp)

def read_acceptor(fname):
    """
    Reads ppm model (e.g. out of 100), and converts to pfm (0.0 - 1.0)
    then creates a dictionary mapping base to a list of frequencies at each
    position. This can then be accessed by the base being scored, and it's 
    position in the motif.
    """
    acc_dict = {}
    f = open(fname)
    for line in f:
        # there are some empty lines at the end of the file??
        # should sort that out in the model script...
        if line != "\n":
            line = line.split()
            base = line[0]
            acc_dict[base] = []
            acc_dict["N"] = []
            for i in range(1,len(line)):
                # add pseudocounts to model, to avoid multiply score by 0
                acc_dict[base].append((int(line[i]) + 0.01) / 100.00)
                acc_dict["N"].append(0.1)
    return acc_dict
    
def read_donor(fname):
    """
    Reads ppm model (e.g. out of 100), and converts to pfm (0.0 - 1.0)
    then creates a dictionary mapping base to a list of frequencies at each
    position. This can then be accessed by the base being scored, and it's 
    position in the motif.
    """
    don_dict = {}
    f = open(fname)
    for line in f:
        # there are some empty lines at the end of the file??
        # should sort that out in the model script...
        if line != "\n":
            line = line.split()
            base = line[0]
            don_dict[base] = []
            don_dict["N"] = []
            for i in range(1,len(line)):
                # add pseudocounts to model, to avoid multiply score by 0
                don_dict[base].append((int(line[i]) + 0.01) / 100.00)
                don_dict["N"].append(0.1)
    return don_dict
    
def score(sequence, donor_model, acceptor_model):
    """
    When given a sequence of the correct length (model to use is determined
    by sequence length), applies a multiplicative scoring of the sequence against
    that model, and returns the score.
    A perfect score would be a model where each position is a completely
    invariant base (i.e. frequency of 1) and every position in the sequence
    matches that model. This would score 1. Anything else is likely to be 
    considerably lower. e.g. perfect donor consensus (from S&S) = 0.29387, 
    perfect acceptor consensus = 0.0001 (because it is more tolerant of change)
    """

    # load the correct model for the input length
    if len(sequence) == 9:
        model = donor_model
    elif len(sequence) == 15:
        model = acceptor_model
    else:
        print "ERROR: Incorrect sequence length"
    
    # apply scoring - multiply scores together
    # pseudocounts ensure that 0 frequency bases are heavily penalised but don't
    # give an automatic 0 score
    score = None
    for i in range(len(sequence)):
        base = sequence[i]
        if score == None:
            score = model[base][i]
        else:
            score = score * model[base][i]
    # consensus acceptor scores a lot lower than consensus donor, so multiply
    # acceptor scores by 300 to normalise.
    if len(sequence) == 15:
	    return score * 300        
    return score
    
def load_models():
    """
    loads the donor and acceptor pwm models
    """
    
    donor_model = read_donor("../shapiro_senepathy_comparison/ppm.canonical.donors")
    acceptor_model = read_acceptor("../shapiro_senepathy_comparison/ppm.canonical.acceptors")
    
    return donor_model, acceptor_model
    
def run():
    """
    Input file should be in the same format as for MESpy, but with the acceptor
    length changed to fit.
    
    For each line of the file, calculate the wt and variant scores, then print
    along with the ID and change in score (absolute and %). Include a header.
    """
    # use sys to get target file from args
    import sys
    
    # load the models
    donor_model, acceptor_model = load_models()
    
    # print a header
    print "ID\tTYPE\tWT\tVAR"    
    
    f = open(sys.argv[1])
    for line in f:
        line = line.rstrip().split()
        var_id = line[1]
        if line[5] == "-":
            wt = reversecomplement(line[2])
            var = reversecomplement(line[3])
        else:
            wt = line[2]
            var = line[3]
        wt_score = score(wt, donor_model, acceptor_model)
        var_score = score(var, donor_model, acceptor_model)
        site_type = line[4]
        print "%s\t%s\t%.9f\t%.9f" % (var_id, site_type, wt_score, var_score)

if __name__ == "__main__":
    
    run()
    
    """
    test_don = "AGTCAATGC"
    test_acc = "AGTCCGACTAGCAAT"
    con_don = "CAGGTAAGT"
    con_acc = "TTTTTTTTTTTCAGG"
    rev_acc = "CCTGAAAAAAAAAAA"
    
    print score(test_don, donor_model, acceptor_model)
    print "--------------"
    print score(test_acc, donor_model, acceptor_model)
    print "--------------"
    print score(con_don, donor_model, acceptor_model)
    print "--------------"
    print score(con_acc, donor_model, acceptor_model)
    """
