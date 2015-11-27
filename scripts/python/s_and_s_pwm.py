import sys

def makepwm(fname):
    """
    Reads a file containing a list of sequences (all same length, no ID fields) and generates a position weight matrix for the sequences.
    This is, at this stage, a list of dicts (one per position) containing the counts of each base at that position across all samples.
    """
    
    f = open(fname, "r")
    nt_dict = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    sequence_list = []  
    
    # use to process the first line on its own  
    first_seq = True
    
    for line in f:
        line = line.strip()
        if "N" in line:
            pass
        else:
            if first_seq:
                # as this is the first sequence, add the dictionary to each position in the list
                i = 0
                for base in line:
                    # add an empty dict to each position in the list
                    sequence_list.append({"A": 0, "C": 0, "G": 0, "T": 0, "N": 0})
                    sequence_list[i][base] += 1
                    i += 1
                # mark as past first line
                first_seq = False
            else:
                i = 0
                for base in line:
                    # just update the list
                    sequence_list[i][base] += 1
                    i += 1
    return sequence_list
        
def percentpwm(rawpwm):
    """
    Takes the dict-based pwm output of makepwm() and converts the raw counts into percentages.
    As there is no dynamic modification needed, the results are returned as a list of tuples,
    ordered (A, C, G, T).
    """
    
    pcpwm = []
    for line in rawpwm:
        # NOTE: This is ignoring "N" bases completely.
        total = line["A"] + line["C"] + line["G"] + line["T"]
        a = (round(line["A"]/float(total),2))*100
        c = (round(line["C"]/float(total),2))*100
        g = (round(line["G"]/float(total),2))*100
        t = (round(line["T"]/float(total),2))*100
        pcpwm.append((a, c, g, t))
    return pcpwm
    
    
def printpwm(pcpwm):
    """
    writes the percentage pwm to the stdout, formatting each number to 3 digits, including 
    leading 0s (000-100), to allow each position to be neatly aligned.
    """
    
    labels = ["A", "C", "G", "T"]
    
    for i in range(4):    
        print labels[i],
        for line in pcpwm:
            print "%03d" % line[i],
        print ""
    print ""
    
def printforlogo(pcpwm):
    """
    Sequence logo generators don't handle percentage inputs, they work from passing in the raw sequences.
    This analysis has far too many sequences for this, so this function generates fake sequences
    that will give the same proportions of bases at each position.
    """
    logo = []
    for col in pcpwm:
        logocol = []
        for i in range(100):
            if i < col[0]:
                logocol.append("A")
            elif i < col[0]+col[1]:
                logocol.append("C")
            elif i < col[0]+col[1]+col[2]:
                logocol.append("G")
            else:
                logocol.append("T")
        logo.append(logocol)
        
    for j in range(100):
        print ">reconstructed_frequency_sequence_%d" % j
        rfseq = []
        for col in logo:            
            rfseq.append(col[j])
        print "".join(rfseq)
        
if __name__ == "__main__":
    # these two arrays are representations of the frequencies from the original
    # shapiro and senapathy models, so they can be recreated as sequence logos.
    sandsdonor = [(32, 37, 19, 12),
                  (58, 13, 15, 15),
                  (10, 4, 78, 8),
                  (0, 0, 100, 0),
                  (0, 0, 0, 100),
                  (57, 2, 39, 2),
                  (71, 8, 12, 9),
                  (5, 6, 87, 5),
                  (16, 15, 22, 47)]

    sandsacceptor = [(9, 31, 15, 45),
                     (9, 33, 13, 45),
                     (7, 31, 11, 51),
                     (7, 35, 7, 51),
                     (10, 35, 7, 47),
                     (10, 35, 11, 44),
                     (7, 43, 7, 42),
                     (9, 41, 8, 42),
                     (6, 39, 6, 48),
                     (6, 40, 8, 46),
                     (23, 29, 23, 27),
                     (3, 74, 1, 22),
                     (100, 0, 0, 0),
                     (0, 0, 100, 0),
                     (28, 12, 49, 10)]
    # if recreating s and s models, uncomment this and comment out the remainder.                 
    #printforlogo(sandsacceptor)   
                  
    fname = sys.argv[1]
    rawpwm = makepwm(fname)
    
    pcpwm = percentpwm(rawpwm)
    if len(sys.argv) > 2:
        if sys.argv[2] == "--logo":
            printforlogo(pcpwm)
    else:
        printpwm(pcpwm)
    
