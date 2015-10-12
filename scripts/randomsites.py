from random import choice

def gen_sequences(LENGTH, REPEATS):
    """
    generates and prints REPEATS x random DNA sequences of length LENGTH.
    Uses unambiguous bases, no N or other uncertainties.
    
    prints pairs of tab-separated sequences to stdout
    """
    for i in range(REPEATS):
        # clear sequences
        sequence1 = []
        sequence2 = []
        for j in range(LENGTH):
            # generate random bases
            sequence1.append(choice(['A','C','G','T']))
            sequence2.append(choice(['A','C','G','T']))
        # write to stdout
        print ''.join(sequence1),"\t",''.join(sequence2)
    
if __name__ == "__main__":
    gen_sequences(50, 100)
