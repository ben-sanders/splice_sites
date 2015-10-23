import sys

def makepwm(fname):
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
    # quick and dirty
    
    labels = ["A", "C", "G", "T"]
    
    for i in range(4):    
        print labels[i],
        for line in pcpwm:
            print "%03d" % line[i],
        print ""
    print ""

if __name__ == "__main__":
    fname = sys.argv[1]
    rawpwm = makepwm(fname)
    
    pcpwm = percentpwm(rawpwm)
    printpwm(pcpwm)
