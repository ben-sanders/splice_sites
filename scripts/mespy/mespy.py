###############################################################################
#                                                                             #
"""                                                                           
============
mespy - MaxEntScanPy
============

MaxEntScan 3' and 5' (Burge & Yeo) rewritten in Python and combined into 
a single method that scans an input sequence and checks every possible site
to find the highest scoring.

To do: add a threshold cutoff for showing sites (score or number of sites)

-------------------------------------------------------------------------------

Ben Sanders 2015-05-01

""" 
#                                                                             #
###############################################################################

from math import log
from operator import itemgetter
from string import maketrans
import argparse
import sys

def log2(value):
    """actually seems to return the log(value)/log(2), rather than log2(value)"""
    return log(value) / log(2)

def getrest(seq, mode):
    """
    extracts a subset of characters from the given string.
    converts to uppercase and returns
    
    in 5' mode gives - string[0,1,2]3,4[5,6,7,8]
    in 3' mode gives - string[0:18],19,20[21,22,23]
    """
    if mode == 5:
        newstring = seq[0:3] + seq[5:9]
        return newstring.upper()
    elif mode == 3:
        # in score3.pl there is no cast to UPPERCASE. helpful to have it set.
        newstring = seq[0:18] + seq[20:23]
        return newstring.upper()
    else:
        # if the mode is incorrect, fail with a ValueError (or maybe TypeError? - does it really matter?)
        raise ValueError
        
def scoreconsensus(seq, mode):
    """
    sets up 3 scoring matrices.
    Then does something.
    """
    
    # values for matrices are taken from the original Perl script
    # bgd is the same for both, then 5' and 3' scripts have different cons1 and cons2
    bgd   = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
    
    if mode == 5:
        # 5' specific matrices
        cons1 = {"A": 0.004, "C": 0.0032, "G": 0.9896, "T": 0.0032}
        cons2 = {"A": 0.0034, "C": 0.0039, "G": 0.0042, "T": 0.9884}
        # $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]});
        score = (cons1[seq[3]] * cons2[seq[4]]) / (bgd[seq[3]] * bgd[seq[4]])
        return score
        
    elif mode == 3:
        # 3' specfic matrices
        cons1 = {"A": 0.9903, "C": 0.0032, "G": 0.0034, "T": 0.0030}
        cons2 = {"A": 0.0027, "C": 0.0037, "G": 0.9905, "T": 0.0030}
        # $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]});
        score = (cons1[seq[18]] * cons2[seq[19]]) / (bgd[seq[18]] * bgd[seq[19]])
        return score
        
    else:
        # if the mode is incorrect, fail with a ValueError (or maybe TypeError? - does it really matter?)
        raise ValueError

def makescorematrix(fname):
    """
    used in 5' ss analysis.
    
    Makes a matrix used for scoring something.
    I'm still figuring out exactly what.
    And also how.
    
    Uses a dictionary in the original, copied here. But should try
    to rewrite using a list to make it a bit more efficient
    """
    try:
        SCOREF = open(fname, "r")
    except:
        print "Can't open %s!" % fname
    matrix = {}
    i = 0
    for line in SCOREF:  
        # remove newlines and whitespace
        line = line.strip()
        line.replace(" ", "")
        # add to dict with i as key, cast to float as I'm pretty sure we want a number!
        matrix[i] = float(line)
        i += 1
    SCOREF.close()
    return matrix
    
def makesequencematrix(fname):
    """
    used in 5' ss analysis
    
    similar to makescorematrix, but works from a list of (every possible?) 7bp
    sequence. creates a dict (this time it actually needs a dict) with the 
    sequence as key and i as the value.
    """
    try:
        SCOREF = open(fname, "r")
    except:
        print "Can't open %s!" % fname
    matrix = {}
    i = 0
    for line in SCOREF:  
        # remove newlines and whitespace
        line = line.strip()
        line.replace(" ", "")
        # add to dict with i as key, cast to float as I'm pretty sure we want a number!
        matrix[line] = i
        i += 1
    SCOREF.close()
    return matrix

def makemaxentscores():
    """
    used in 3' ss analysis
    
    generates a list of models based on the files below.
    not sure what these models do or how, but each model is a dict
    of index : sequence.
    """
    modeldir = sys.path[0]+"/splicemodels/";
    modellist = ['me2x3acc1', 'me2x3acc2', 'me2x3acc3',
                 'me2x3acc4', 'me2x3acc5', 'me2x3acc6',
                 'me2x3acc7', 'me2x3acc8', 'me2x3acc9']
    metables = []
    for model in modellist:
        # try to open the model file
        try:
            SCOREF = open(modeldir + model)
        except:
            print "Can't open %s!" % model
        
        # each model gets a dictionary, key i and value line
        # these are stored in the metables array
        n = 0
        modeldict = {}
        for line in SCOREF:
            line = line.strip()
            line = line.replace(" ", "")
            modeldict[n] = line
            n += 1
            
        # add dict for this model to the list, then close the file
        metables.append(modeldict)
        SCOREF.close()
    return metables

def hashseq(sequence):
    """
    used in 3' ss analysis.
    
    wow! documentation from the original!
    
    returns hash of sequence in base 4
    hashseq('CAGAAGT') returns 4619
    """    
    sequence = sequence.translate(maketrans("ACGT", "0123"))   
    seqsum = 0
    length = len(sequence)
    four = [1, 4, 16, 64, 256, 1024, 4096, 16384]
    for i in range(length):
        # $sum+= $seqa[$i] * $four[$len - $i -1];
        # sequence.translate used chars not int, so need to cast the sequence
        seqsum += int(sequence[i]) * four[length - i -1]
    return seqsum

def maxentscore(sequence, tables):
    """
    Used in 3' ss analysis.
    
    generates a score by looking up in a table.
    I don't actually know what the tables represent, or how this works,
    but it's just adapted from the Perl script
    """
    
    sc = []
    # I have no idea what these various tables and lookups actually mean.
    sc.append(float(tables[0][hashseq(sequence[0:0+7])]))
    sc.append(float(tables[1][hashseq(sequence[7:7+7])]))
    sc.append(float(tables[2][hashseq(sequence[14:14+7])]))
    sc.append(float(tables[3][hashseq(sequence[4:4+7])]))
    sc.append(float(tables[4][hashseq(sequence[11:11+7])]))
    sc.append(float(tables[5][hashseq(sequence[4:4+3])]))
    sc.append(float(tables[6][hashseq(sequence[7:7+4])]))
    sc.append(float(tables[7][hashseq(sequence[11:11+3])]))
    sc.append(float(tables[8][hashseq(sequence[14:14+4])]))
    score = sc[0] * sc[1] * sc[2] * sc[3] * sc[4] / (sc[5] * sc[6] * sc[7] * sc[8])
    return score

def openfile(fname):
    """
    opens the given file, and concatenates all lines into one single sequence.
    if there is a fasta descriptor in the first line, it is skipped.
    if there is a fasta descriptor anywhere else in the file input is stopped,
    and only the sequence up to that point is used.
    
    ignoring fasta descriptor, and using file name as the identifier. Why? 
    just because it's simpler.
    """
     
    # open the file
    f = open(fname, "r")
    # merge the file sequence into one string
    # lc counts lins - if a ">" appears on a line other than 0, stop
    lc = 0
    sequence = ""
    for line in f:
        if line.startswith(">"):
            if lc != 0:
                f.close()
                break
            else:
                lc += 1
                pass
        else:
            lc += 1
            # remove newline
            line = line.strip()
            # add line to sequence so far
            sequence = sequence + line
    # ensure the sequence is uppercase
    sequence = sequence.upper()
    return sequence
  
def scorefive(me2x5, seq, sequence):
    """scores sequence with a 5' splice site model"""
    #,&log2(&scoreconsensus($str)*$me2x5{$seq{&getrest($str)}})
    score = log2(scoreconsensus(sequence, 5) * me2x5[seq[getrest(sequence, 5)]])
    return score

def scorethree(sequence, tables):
    """scores sequence with a 3' splice site model"""
    #log2(scoreconsensus($str)*maxentscore(getrest($str),\@metables))
    score = log2(scoreconsensus(sequence.upper(), 3) * maxentscore(getrest(sequence, 3), tables))
    return score
    
def window(sequence, mode="both"):
    """
    run along sequence and split into windows of the correct sizes
    for 5' 9nt
    for 3' 23nt
    
    in each window, call the appropriate score function and record result.
    """
    
    # create the score matrices here, or they'll be recreated with each score call
    me2x5 = makescorematrix(sys.path[0]+'/splicemodels/me2x5')
    seq = makesequencematrix(sys.path[0]+'/splicemodels/splice5sequences')
    threetables = makemaxentscores()
    
    scores = []
    
    if mode == 5 or mode == "both":
        # 5' windowing, size 9nt
        windowsize = 9
        for i in range(len(sequence) - windowsize + 1):
            window = sequence[i:i+windowsize]
            # format for display, lowercase intronic, uppercase exonic
            printwindow = window[:3].upper() + window[3:].lower()            
            scores.append((5, i+1, printwindow, scorefive(me2x5, seq, window)))
    if mode == 3 or mode == "both":
        # 3' windowing, size 23nt
        windowsize = 23
        for i in range(len(sequence) - windowsize + 1):
            window = sequence[i:i+windowsize]
            # format for display, lowercase intronic, uppercase exonic
            printwindow = window[:20].lower() + window[20:].upper()      
            scores.append((3, i+1, printwindow, scorethree(window, threetables)))
    
    return scores
        
def run(fname, threshold, top, verbose, exclusive):
    """
    wrapper to run the whole analysis
    """
    # load the sequence    
    sequence = openfile(fname)
    # get the scores  
    scores = window(sequence)
    # sort the scores
    scores = sorted(scores, key=itemgetter(3), reverse=True)
    
    # check if any output filters are specified                
    # showing a quick summary if verbose mode is True
    # borrowing vcf style line indicators: header = ##, column names = #, data blank
    if verbose == True:
        print "##-------------------------------------------------------------"
        print "## mespy - MaxEntScanPy results for %s\n##" % fname
        print "## %d potential splice sites analysed." % len(scores)
        if exclusive != False:
            print "## Showing %d' sites only." % exclusive
        if top != False:
            print "## Showing top %d results only." % top
        if threshold != False:
            print "## Showing scores above %.2f only." % threshold
        print "##-------------------------------------------------------------\n##"
    
    # quick sanity check on exclusive site filter
    # warning will display even in non-verbose mode
    if exclusive not in [3, 5, False]:
        print "## WARNING: exclusive splice site filter not recognised. Showing both 3' and 5' sites"
        exclusive = False
    
    # always show column names
    print "#Type\tStart\tSequence\tScore"
    
    # filter by exclusive site type
    if exclusive != False:
        outscores = []
        for score in scores:
            if score[0] == exclusive:
                outscores.append(score)
    else:
        outscores = scores[:]
    
    # filter by <n> top scoring results
    if top != False:    
        if top <= len(scores): # just in case someone enters too high a cutoff
            outscores = outscores[:top]
        else:
            outscores = outscores[:]
    else:
        outscores = outscores[:]
        
    # filter by score - already sorted, so print as it goes
    if threshold != False:
        for score in outscores:
            if score[3] >= threshold:
                print "%d\t%d\t%s\t%.3f" % (score[0], score[1], score[2], score[3])
                
    # if no score filter, print everything in outscores
    else:
        for score in outscores:
            print "%d\t%d\t%s\t%.3f" % (score[0], score[1], score[2], score[3])
   
if __name__ == "__main__":
    # configure the argparser
    parser = argparse.ArgumentParser(description="""mespy - MaxEntScanPy. Scores potential splice sites in a sequence file""")
    parser.add_argument("fname", metavar="<f>", type=str, help="Name of file to process")
    parser.add_argument("-x", "--exclusive", type=int, default=False, help="show one type of site exclusively. False = both")
    parser.add_argument("-s", "--score", type=float, default=False, help="Set score cutoff.")
    parser.add_argument("-t", "--top", type=int, default=False, help="Display top n sites only.")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose mode.")
    
    args = parser.parse_args()
    
    # go!
    run(args.fname, args.score, args.top, args.verbose, args.exclusive)
