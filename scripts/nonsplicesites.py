# go through the hg19.fa and fine all defined dinucleotide positions in 
# bed format

def run(fname, target):
    # open genome fasta for reading
    f = open(fname, "r")
    
    # track the current chromosome and position - for bed output
    currentchrom = None
    lastbase = None
    currentbase = None
    currentpos = 0
    counter = 1
    
    for line in f:
        # remove trailing newlines
        line = line.strip()
        if line.startswith(">"):
            # delete the > later
            currentchrom = line[1::]
            # reset base trackers at start of new chromosome
            lastbase = None
            currentbase = None
            currentpos = 0
        else:
            for i in range(len(line)):
                # track the current, and previous bases, then concatenate and
                # compare against the target
                if currentbase == None:
                    # if this is the first base, can't check a dinucleotide
                    # just record current base and move to next.
                    currentbase = line[i]
                else:
                    # move the last base, and get the current one.
                    # then compare against the target
                    lastbase = currentbase
                    currentbase = line[i]
                    window = lastbase + currentbase
                    if window == target:
                        # output in bed format, 
                        if target == "AG":
                            print "%s\t%d\t%d\tRAND_%06d_acc\t1\t+" % (currentchrom, 
                                                                       currentpos+1, 
                                                                       currentpos+1,
                                                                       counter)
                        else:
                            print "%s\t%d\t%d\tRAND_%06d_don\t0\t+" % (currentchrom, 
                                                                       currentpos-1, 
                                                                       currentpos-1,
                                                                       counter)
                        counter += 1
                # increment the position (do at end to retain 0-indexing)
                # doing at beginning would mean starting at 1, which is not
                # right for bedtools.
                currentpos += 1

if __name__ == "__main__":
    import sys
    # could test inputs, but will be called from a script so that's controlled
    genomefile = sys.argv[1]
    dinucleotide = sys.argv[2]
    run(genomefile, dinucleotide)
