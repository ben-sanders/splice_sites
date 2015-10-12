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
