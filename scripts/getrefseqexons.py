"""
---------------------------
--<< getrefseqexons.py >>--
---------------------------
Ben Sanders
2015-08-19

inputs:
    txt file output from UCSC table browser
        group: Genes and Gene Predictions
        track: RefSeq Genes
        table: refGene

outputs:
    to stdout, one transcript per gene (with most exons) in BED format
        chromosome
        start
        end (should be same as start)
        <gene>_int<intron_num>_[don/acc]
        type (0 = donor, 1 = acceptor)
        strand [+/-]

"""

class Gene(object):
    """
    Store the information from a gene in the RefSeq genes list.    
    Need to split this data up into a list of exon start/finish coordinates.    
    """
    
    def __init__(self, name):
        self.transcriptlist = {}
        self.name = name
        
    def addtranscript(self, transcript):
        """
        Wrapper to add a transcript to the transcriptlist dict
        """
        self.transcriptlist[transcript.getid()] = transcript
    
    def getname(self):
        """
        Returns the gene name, which actually corresponds to the name2 field
        from the UCSC data. Should be the HGNC official symbol.
        """
        return self.name
        
    def getlongest(self):
        """
        Returns the transcript with the most exons.
        This is not necessarily the most common transcript, but it has the
        most splice sites in it, which is what I'm after.
        """
        # variables to track the longest transcript
        self.mostexons = 0
        self.mosttranscript = None
        # each transcript has an exonCount variable as part of the UCSC data
        # getexoncount() just returns this value
        for transcript in self.transcriptlist.values():   
            if transcript.getexoncount() >= self.mostexons:
                self.mostexons = transcript.getexoncount()
                self.mosttranscript = transcript
        return self.mosttranscript
    
    def getsplicesites(self):
        """
        gets the longest transcript to print out a list of it's splice sites
        in BED format.
        """
        # get the start and end coordinates for the transcript with the most
        # exons. Not necessarily the most frequent splice sites, but they
        # are biologically functional (hopefully...)
        self.longesttranscript = self.getlongest()
        self.longesttranscript.getdata()
            
    def __str__(self):
        print self.name
        for transcript in self.transcriptlist.keys():
            print "\t",transcript
        return "%s: %d" % (self.getlongest().getid(), self.getlongest().getexoncount())
        
class Transcript(object):
    """
    store an individual transcript, so we can keep them together within a gene.
    """
    def __init__(self, line):
        line = line.rstrip().split()
        self.bin        = line[0]
        self.name       = line[1]
        self.chrom      = line[2]
        self.strand     = line[3]
        self.txStart    = int(line[4])
        self.txEnd      = int(line[5])
        self.cdsStart   = int(line[6])
        self.cdsEnd     = int(line[7])
        self.exonCount  = int(line[8])
        # exon start and end coordinates are comma-separated lists
        # there is a trailing comma, which leads to a blank entry at the end
        # of the post-split list, so remove the last entry
        self.exonStarts = line[9].split(",")[:-1]
        self.exonEnds   = line[10].split(",")[:-1]
        self.score      = int(line[11])
        # name2 is the actual gene symbol
        self.name2      = line[12]
        self.cdsStartStat = line[13]
        self.cdsEndStat = line[14]
        self.exonFrames = line[15]
        
    def getgene(self):
        """
        return the gene name, so the transcript can be added to the Gene 
        record.
        """
        return self.name2
        
    def getid(self):
        return self.name
        
    def getexoncount(self):
        return self.exonCount
        
    def getexons(self):
        return (self.exonStarts, self.exonEnds)
        
    def getstrand(self):
        return self.strand
        
    def getchrom(self):        
        return self.chrom
        
    def getdata(self):        
        """
        break the longest trascript into a list of splice site positions.
        split the exonStarts and exonEnds lists from the UCSC data, and order
        according to strand direction.
        Remember to always account for direction! splicing is relative to the
        transcribed RNA, so need to reverse (and complement?) sequences from 
        negative strand
        
        # REMEMBER: structure is
        #   5'---exon---donor|---intron---bs--ppt--|acceptor---intron---3'
        # so on +ve self.starts are acceptors.
        
        """  
        # do a quick sanity check on the positions
        # data is given relative to + strand, so all values should increase,
        # regardless of the actual gene orientation.
        # e.g. for - strand genes, exonStarts are actually exonEnds
        # This is sorted later.
        assert int(self.exonStarts[0]) < int(self.exonStarts[-1]), "?reverse strands"+self.name2
        assert int(self.exonStarts[0]) < int(self.exonEnds[0]), "?reverse strands"+self.name2
        
        # if the strand is +, the orientation is ok as it is
        if self.strand == "+":
            # starts = donors, ends = acceptors
            # ignore the first start and last end, as these aren't splice sites
            self.acceptors = self.exonStarts[1:]
            self.donors = self.exonEnds[:-1]
        # negative strand these should be reversed, as positions are given
        # relative to + strand
        else:
            # gene is on negative strand, but numbered by positive strand.
            # reverse the lists to correct the orientation
            self.exonStarts.reverse()
            self.exonEnds.reverse()
            # as the gene is on the reverse strand, donors and acceptors are
            # reversed - i.e. exonEnds are acceptors, exonStarts are donors
            # as before, ignore the non-splice site exon boundaries.
            self.acceptors = self.exonEnds[1:]
            self.donors = self.exonStarts[:-1]
            
        # let's have some more sanity checks
        # order should now be acceptors[0], donors[1], acceptors[1],
        # because the start of exon 1 is the transcription start, not a splice
        # site. Also the end of the last exon is transcription stop, not ss.
        # orientation matches
        assert len(self.donors) == len(self.acceptors), "?number of sites"
        if self.strand == "+":            
            if len(self.donors) > 1 and len(self.acceptors) > 1:
                assert int(self.acceptors[0]) < int(self.acceptors[-1]), "?correct strand"+self.name2
                # first donor should be before first acceptor (see diagram above)
                assert int(self.donors[0]) < int(self.acceptors[0]), "?site ordering"+self.name2
        else:
            if len(self.donors) > 1 and len(self.acceptors) > 1:   
                assert int(self.acceptors[0]) > int(self.acceptors[-1]), "?correct strand"+self.name2
                assert int(self.donors[0]) > int(self.acceptors[0]), "?site ordering"+self.name2
            
        # Print in BED format, to get sequence with bedtools.
        # write to stdout, so can redirect if wanted.
        for i in range(len(self.donors)):
            # print donor as a single base
            print "%s\t%s\t%s\t%s_int%d_don_%s\t0" % (self.chrom, 
                                                         self.donors[i],
                                                         self.donors[i],
                                                         self.name2,
                                                         i+1,
                                                         self.strand)
            # now acceptor - set score to 1
            print "%s\t%s\t%s\t%s_int%d_acc_%s\t1" % (self.chrom, 
                                                         self.acceptors[i],
                                                         self.acceptors[i],
                                                         self.name2,
                                                         i+1,
                                                         self.strand)
            
    def __str__(self):
        return "%s\t%s\t%s\t%s" % (self.chrom, self.name2, self.name, self.exonEnds)

class GeneList(object):
    """
    store a dict of gene names, link to a list (or another dict?) of each 
    associated transcript
    """
    def __init__(self):
        self.genelist = {}
        
    def addtranscript(self, transcript):
        self.allowedchroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                         "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                         "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                         "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        # if the gene already exists, add the transcript to its dict
        # if not, create the gene and add the transcript
        
        # only add a transcript (or create a gene) if there are multiple exons
        # otherwise there are no splice sites, and what use is that!?
        # also don't want unecessary creation of objects - i.e. useless
        # Transcripts and empty Genes.
        if not transcript.getexoncount() <= 1:
            if transcript.getchrom() in self.allowedchroms:
                if transcript.getgene() not in self.genelist.keys():
                    self.genelist[transcript.getgene()] = Gene(transcript.getgene())
                self.genelist[transcript.getgene()].addtranscript(transcript)
        
    def getsplicesites(self):
        """
        For each gene in the list, runs the Gene.getsplicesites() method to 
        return a list of splice site locations within that gene.
        """
        # for now this just prints on each call, but in the future
        # should perhaps store in a list?
        # or stick with print so output can be redirected.
        for gene in self.genelist.keys():
            self.genelist[gene].getsplicesites()
        
    def getgene(self, name):
        """
        Returns a given Gene object, or returns an error.
        Could raise a KeyError, but that might break things?
        """
        return self.genelist.get(name, "ERROR: Gene not found")
        
    def __str__(self):
        #    print gene  
        #    print self.genelist[gene]
        #    for transcript in self.genelist[gene].keys():
        #        print "\t", transcript
        for gene in self.genelist.keys():
            print self.genelist[gene]
        return "----------------------"
    
def run(fname):
    # pretty much all the work is done within GeneList, Gene, and Transcript 
    # objects
    genelist = GeneList()
    
    # each line is a transcript, parsed and filtered by the Gene to which the
    # GeneList assigns it.
    f = open(fname, "r")
    for line in f:
        if not line.startswith("#"):
            # cast line to a Transcript and send to GeneList.
            genelist.addtranscript(Transcript(line))
    
    # works out and prints splice sites
    # might make it to return a list, or just rejig output for redirecting
    genelist.getsplicesites()

if __name__ == "__main__":
    run("refgenes.txt")
