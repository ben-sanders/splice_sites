# reads a list of FASTA sequences
# if sequence is on reverse strand (will be in identiifier, reverse complement
# read sequence by sequence, and print to stdout.

import sys  # to write to stderr

def read(fname):
    f = open(fname, "r")
    to_rev_comp = False
    linecount = 0
    revcompcount = 0
    
    for line in f:
        linecount += 1
        if to_rev_comp:
            # reverse complement + print this line
            print reversecompliment(line.rstrip())
            revcompcount += 1
            # reset the revcomp marker
            to_rev_comp = False
        # identifiers preceding lines to revcomp end in "-"
        elif line.startswith(">") and line.endswith("-\n"): 
            # revcomp the next line (containing sequence)...
            to_rev_comp = True
            # print identifier	
            print line.rstrip()
        else:
            print line.rstrip()
    
    # summary should be pritned to stderr, so stdout can be redirected to file without it.
    name = "-------------------\n" + str(fname) + "\n-------------------\n"
    sys.stderr.write(name)
    tl = "total lines: " + str(linecount) + "\n"
    sys.stderr.write(tl)
    ts = "total sequences: " + str(linecount / 2) + "\n"
    sys.stderr.write(ts)
    rc = "rev comped: " + str(revcompcount) + "\n"
    sys.stderr.write(rc)     
    
def reversecompliment(sequence):
    """
    Creates the reverse compliment of a DNA sequence
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
    
if __name__ == "__main__":
    read(sys.argv[1])


