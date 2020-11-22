#!/usr/bin/env python3
# Name: Matthew Kozubov(mkozubov)
# Group Members: Liana Beld

'''
Desktop/python scripts/scripting/Todd_intron_script/132_Todd_introns.txt
'''

import sys

'''
My bioinformatics tool kit.
* Fasta reader can read in fasta files
* NucParams can give codon composition info
* Protein params, gives information on a peptide sequence.
* ORFfinder finds possible, non-intronic, genes
* UniqueSubstring can generate unique and essential substrings
  of a list of strings
* tRNA class can be modified to hold more info and generate a better
  tRNA model.

'''

class ORFfinder:
    '''
    This class finds possible genes within a given raw sequence
    on both the bottom and top strand, checking +/-1,2,3 frames.
    Returns a list of tuples. Indexing is 1-len(seq)
    
    Has the optional parameters: 
    * minSize (genes need to be at least this big to be considered)
    * largestORF (only the largest start-stop pair will be considered)
    * startCodons (a mutable set of possible starts)
    * stopCodons (a mutable set of possible stops)
    
    initialization:
    test = ORFfinder('sequence')

    usage:
    test.getOrfs()

    return:
    [(frame, start, stop, length)...]

    '''
    complementDic = {'A':'T',
                     'T':'A',
                     'C':'G',
                     'G':'C'}

    def __init__(self, seq, minSize = 100, largestORF = True, startCodons = 'ATG', stopCodons = ['TAG', 'TAA', 'TGA']):
        '''constructor: initializes passed optional parameters for use in the _findOrfs method'''
        self.seq = seq.upper()
        self.minSize = minSize
        self.largestORF = largestORF
        self.startCodons = startCodons
        self.stopCodons = stopCodons

        self.orfList = []


    def _reverseComplement(self):
        '''Returns the reverse complement in all-caps'''
        
        #found from www.educative.io, how to reverse a string
        seq = self.seq
        rseq = ''.join(reversed(seq.upper()))
        rcseq = ''
        for nuc in rseq:
            rcseq += self.complementDic[nuc]
        return rcseq

    def _findOrfs(self, seq, reverse = False):
        '''Populates the orfList with tuples of (frame, start position, stop position, length). Checks left, right, and dangling edge cases.'''
        if reverse:
            seq = self._reverseComplement()
        
        for frame in range(0,3):
            startCodonList = []
            stopCodonCounter = 0
            #handles middle cases
            for index in range(frame, len(seq), 3):
                #codon was filled with fragments, messing up my in operation, thus the if statement
                codon = seq[index:index+3] if len(seq[index:index+3]) == 3 else 'not'

                if codon in self.startCodons:#ask if codon is a start
                    startCodonList.append(index)
                elif codon in self.stopCodons:#ask if codon is a stop
                    stopCodonCounter += 1
                    if self.largestORF and startCodonList:#asks if we want the largest orf
                        if index+3-startCodonList[0] > self.minSize:#asks if the gene too append is larger than minSize
                            if reverse:#asks if this is bottom strand
                                self.orfList.append((-(frame+1), (len(seq)-index-3)+1, (len(seq)-startCodonList[0]-1)+1, index+3-startCodonList[0]))
                            else:#asks if this is top strand
                                self.orfList.append((frame+1, 1+startCodonList[0], 1+index+2, index+3-startCodonList[0]))
                    else:#in the case we want all orfs
                        for start in startCodonList:
                            if index+3-start > self.minSize:#asks if the gene too append is larger than minSize
                                if reverse:#asks if this is the bottom strand
                                    self.orfList.append((-(frame+1), 1+len(seq)-index-3, 1+len(seq)-start-1, index+3-start))
                                else:#asks if this is the top strand
                                    self.orfList.append((frame+1, 1+start, 1+index+2, index+3-start))
                    startCodonList = []

            #handles dangling boys
            #asks if we found any start codons/stop codons
            if not startCodonList and not stopCodonCounter:
                if len(seq) > self.minSize:#asks if the gene too append is larger than minSize
                    if reverse:
                        self.orfList.append(( -(frame+1), 1, len(seq), len(seq)))
                    else:
                        self.orfList.append((frame+1, 1, len(seq), len(seq)))

            
            #if we exit out of the loop with startcodons, but no stop
            #handles right edge case
            if startCodonList:#asks if we found start codons, but no stop codons
                for start in startCodonList:
                    if len(seq)-start > self.minSize:#asks if the gene too append is larger than minSize
                        if reverse:#asks if bottom strand
                            self.orfList.append((-(frame+1), 1, len(seq)-start, len(seq)-start))
                        else:#asks if top strand
                            self.orfList.append((frame+1, start+1, len(seq), len(seq)-start))
                startCodonList = []  

            #handles left edge case
            for index in range(frame, len(seq), 3):
                codon = seq[index:index+3] if len(seq[index:index+3]) == 3 else 'not'
                if codon in self.startCodons:#asks if the codon is a start
                    break
                elif codon in self.stopCodons:#asks if the codon is a stop
                    if index+3 > self.minSize:#asks if the gene too append is larger than minSize
                        if reverse:#asks if top strand
                            self.orfList.append( (-(frame+1), len(seq)-index-2, len(seq), index+3) )
                        else:#asks if bottom strand
                            self.orfList.append((frame+1, 1, index+3, index+3))
                        break

    def getOrfs(self):
        '''Returns populated orfList, making use of _findOrfs for top and bottom strand'''
        self._findOrfs(self.seq, reverse = False)#gets top strand
        self._findOrfs(self.seq, reverse = True)#gets bottom strand
        return self.orfList

class NucParams:
    '''
    Given raw DNA or RNA sequences, creates an object that can
    return amino acid composition, codon composition,
    nucleotide composition, and total nuc count.

    It is also mutable, extra sequence can be added via
    the addSequence(str) method.

    initialized:
    param = NucParams(str)
    and or
    param = NucParams()
    param.addSequence(str)

    '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}


    def __init__(self, inString=''):#works
        '''
        constructor: initializes two dictionaries: countDicNuc, countDicaa

        optional: Translates the inString into RNA
                    if it was DNA. Fills the countDicNuc
                    with counts of codons from the
                    passed inString variable. Will chop
                    off hanging nucleotides if they are 
                    not in a multiple of 3
        '''
        self.countDicNuc = {codon:0 for codon in self.rnaCodonTable.keys()}
        self.countDicaa = {aa: 0 for aa in self.rnaCodonTable.values()}
        
        #boolean for nucComp method
        self.rnaQ = True

        if inString:
            self.addSequence(inString)
        else:
            pass

    def addSequence(self, inSeq):#works, heavy lifter
        '''
        Takes an input string of {ATGCNU} characters
        and appends it to the existing countDicNuc 
        dictionary.

        Returns: nothing
        '''
        #checks if RNA or DNA
        if inSeq.find('T') > -1:
            self.rnaQ = False

        inSeq = inSeq.upper().replace(' ', '').replace('T', 'U')

        for start in range(0, len(inSeq), 3):
            triplet = inSeq[start:start+3]
            if len(triplet)%3 != 0:
                break
            #can handle N's
            self.countDicNuc[triplet] = (self.countDicNuc.get(triplet) if self.countDicNuc.get(triplet) else 0) + 1

    def aaComposition(self):#works
        '''Translates the dictionary of nuc counts into aa. Returns the aa dictionary'''

        for nuc, count in self.countDicNuc.items():
            #can handle N's
            if nuc in self.rnaCodonTable.keys():
                aa = self.rnaCodonTable.get(nuc)
                self.countDicaa[aa] += count

        return self.countDicaa

    def codonComposition(self):#works
        """Returns a dictionary of valid codons, ignores N's and illegal characters"""
        return {codon: self.countDicNuc.get(codon) for codon in self.rnaCodonTable.keys()}
    
    def nucComposition(self):#works
        '''Returns a dictionary of valid nucleotides and their counts'''
        if self.rnaQ:
            #is rna
            rnaDic = {'A': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}

            for key, n in self.countDicNuc.items():
                A = n*key.count('A');C = n*key.count('C');G = n*key.count('G')
                U = n*key.count('U');N = n*key.count('N')

                rnaDic['A'] += A; rnaDic['C'] += C; rnaDic['G'] += G
                rnaDic['U'] += U; rnaDic['N'] += N

            return rnaDic

        else:
            #is dna
            dnaDic = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
            
            for key, n in self.countDicNuc.items():
                A = n*key.count('A');C = n*key.count('C');G = n*key.count('G')
                T = n*key.count('U');N = n*key.count('N')

                dnaDic['A'] += A; dnaDic['C'] += C; dnaDic['G'] += G
                dnaDic['T'] += T; dnaDic['N'] += N

            return dnaDic

    def nucCount(self):#works
        '''Returns the number of valid nucleotides in the given sequence as an int'''
        return sum(count for count in self.nucComposition().values())

class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        #yield makes this compatable with a for loop
        yield header,sequence

class ProteinParam:
    '''
    This class can return amino acid counts, theoretical pI, amino acid composition,
    molar and mass extinction coefficients, and molecular weight of an input
    peptide sequence.
    
    initialization:
    test = ProteinParam(sequence)

    usage:
    test.aaCount(); test.pI(); test.aaComposition();
    test.molarExtinction(); test.massExtinction(); 
    test.molecularWeight()
    '''
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        initializes the object, converting the input string to uppercase
        as well as creating a dictionary of aa counts.
        '''
        self.protein = protein.upper()
        self.protDic = {key:(self.protein.count(key)) for key in self.aa2mw.keys()} #works even if there are weird characters!!!!


    def aaCount (self):
        '''returns total number of valid aa as an int'''
        return sum(value for value in self.protDic.values())

    def pI(self):
            '''returns the theoretical pI as a string, precise to 2 decimal places. Works via a binary search. PI = closest num where _charge_ = 0'''

            lowerBoundry = 0
            upperBoundry = 14

            precision = .01

            #range depends on how precise we want to be
            while upperBoundry - lowerBoundry >= precision:

                middlePH = (lowerBoundry + upperBoundry) / 2
                middleCharge = self._charge_(middlePH)

                #checks whether the lower and middle bound have the same sign, positive or negative
                if (middleCharge >= 0):
                    lowerBoundry = middlePH
                  
                else:
                    upperBoundry = middlePH
                    

            return middlePH


    def aaComposition (self) :
        '''returns the dictionary of counts created in the __init__'''
        return self.protDic

    def _isNegative_(self, integer):
        '''returns True if the the passed int is less than 0'''
        return True if integer < 0 else False

    def _charge_(self, pH):
        '''returns the net charge of a protein at a specified pH as a float'''

        aaNterm = self.aaNterm
        aaCterm = self.aaCterm
        aa2chargePos = self.aa2chargePos
        aa2chargeNeg = self.aa2chargeNeg
        protDic = self.protDic

        #math
        positive = sum( (protDic.get(name) if protDic.get(name) else 0) * (10**N/(10**N + 10**pH))\
                for name,N in aa2chargePos.items())

        #math
        negative = sum( (protDic.get(name) if protDic.get(name) else 0) * (10**pH/(10**N + 10**pH))\
                for name,N in aa2chargeNeg.items())

                        #C-terminus                             #N-terminus                         #aa charges
        netCharge = (  ((10**aaNterm)/(10**aaNterm+10**pH)  -  (10**pH)/(10**aaCterm+10**pH))  ) + (positive - negative)


        return netCharge

    def molarExtinction (self):
        '''returns the molar extinction coefficient as an int'''
        num = sum(self.protDic[name]*self.aa2abs280[name] for name in self.aa2abs280.keys()) 
        return num

    def massExtinction (self):
        '''returns the mass extinction coefficient as a float rounded to 2 decimal places'''
        myMW =  self.molecularWeight()
        return round(self.molarExtinction() / myMW, 2) if myMW else 0.0

    def molecularWeight (self):
        '''returns the molecular weight of the protein sequence as a float rounded to 1 decimal place'''
        return round(self.mwH2O+sum(N*(self.aa2mw[name] - self.mwH2O)\
               for name, N in self.protDic.items()), 1)

class UniqueSubstring:
    '''
    Creates an object that, given a list of strings,
    will find the Unique and Essential set of substrings.
    or
    The Ubiquitous and Non-essental set of substrings

    initialization:
    test = UniqueSubstring( list(str) )

    usage:
    test.definingFeatures()
    test.ubiqitiousFeatures()


    '''
    def __init__(self, tRNAList):
        self.tRNAList = tRNAList

    def _strToSet(self, test):#works
        '''Return a set of all contigious combinations given a single string'''
        setOfTest = set()

        for start in range(0,len(test)):
            for stop in range(start+1, len(test)+1):
                setOfTest.add(test[start:stop])

        return setOfTest

    def _listToSet(self, listOftRNAs):#works
        '''Return a list of all contigious combinations given a list of strings'''
        return [self._strToSet(tRNA) for tRNA in listOftRNAs]

    def _unionSet(self, listOfOthers):#works
        '''Return the union of a list of sets'''
        unionSet = set()
        for tRNA in listOfOthers:
            unionSet = unionSet.union(tRNA)
        return unionSet

    def _makeUnique(self, alphaSet):#works
        '''Return a list of UNIQUE sets given the sets of all possible substrings'''
        betaList = []#works
        for n, i in enumerate(alphaSet):
            listOfOthers = alphaSet[:n] + alphaSet[n+1:]
            betaList.append( alphaSet[n].difference( self._unionSet(listOfOthers) ))
        return betaList

    def _makeUbiq(self, alphaSet):
        '''Return a list of UBIQITOUS sets given the sets of all possible substrings'''
        basic = {'A', 'G', 'C', 'U'}
        betaList = []

        for searchNum, searchRNA in enumerate(alphaSet):
            similar = searchRNA
            for currentNum, otherRNA in enumerate(alphaSet):
                if currentNum == searchNum:
                    continue
                similar = similar.intersection(otherRNA)

            betaList.append(similar.difference(basic))

        return betaList

    def _makeEssential(self, betaList):#works better :)
        '''Return a list of ESSENTIAL sets given the list of UNIQUE sets'''
        omegaList = []
        #now need to get ESSENTIAL elements of the set
        for item in betaList:
            alignSet = set()
            for element in item:
                if element[:len(element)-1] in item:
                    alignSet.add( element )
                elif element[1:] in item:
                    alignSet.add( element )

            omegaList.append( item.difference(alignSet) )

        return omegaList

    def _makeEssential_old(self, betaList):#works
        '''Return a list of ESSENTIAL sets given the list of UNIQUE sets'''
        omegaList = []
        #now need to get ESSENTIAL elements of the set
        for item in betaList:
            removeSet = set()
            for element in item:
                for element2 in item:
                    if element is element2:
                        continue
                    elif element in element2:
                        removeSet.add(element2)

            omegaList.append( item.difference(removeSet) )

        return omegaList

    def definingFeatures(self):
        '''Return a list of sets contanining defining features from strings'''
        alphaSet = self._listToSet(self.tRNAList)
        betaList = self._makeUnique(alphaSet)
        omegaList = self._makeEssential(betaList)

        return omegaList

    def ubiqitiousFeatures(self):
        '''Return a list non-essential sets contanining features shared among the strings'''
        alphaSet = self._listToSet(self.tRNAList)
        betaList = self._makeUbiq(alphaSet)
        #omegaList = self._makeEssential(betaList)

        #return omegaList
        return betaList

class tRNA:
    '''
    Given a header and a raw sequence, this class can
    generate the powerset of the sequence, as well as recognize
    what tRNA has been passed.

    initialization:
    test = tRNA(header, seq)

    ###header = 'tRNA|Glu|UGC|Species|moreInfo'
    ###seq = 'GUUCUU...'


    usage:
    sequence = test.seq()
    tRNAname = test.name()
    powerSet = test.powerSet()

    returns:
    sequence = 'GUUCUU...'
    tRNAname = 'Glu'
    powerSet = {'GUUCUU', 'UUCUUG',...}

    '''
    def __init__(self, header, seq):
        '''Initializes the sequence and header'''
        self.seq = seq
        self.head = header
        self.uniqueSet = None

    def sequence(self):
        '''Return the raw sequence'''
        return self.seq

    def header(self):
        '''Returns the header of the tRNA'''
        return self.head

    def name(self):
        '''Return the name of the tRNA'''
        return self.head.split('|')[1]

    def powerSet(self):
        '''Return a set of all contigious combinations of the sequence'''
        setOfTest = set()
        test = self.seq

        for start in range(0, len(test)):
            for stop in range(start+1, len(test)+1):
                setOfTest.add(test[start:stop])

        return setOfTest

    def addUniqueSet(self, genSet):
        '''Given a set, stores that as an attribute in self.uniqueSet'''
        self.uniqueSet = genSet

class FastBoundry:
    '''
    Given a .fa file with one entry, boundries,
    and strand, outputs to stdout the contents
    within the boundries. 5' --> 3' (+) format
    '''
    def __init__(self,  bounds, infile=''):
        self.infile = infile
        self.bounds = bounds

    def printBoundFasta(self):
        ''''''
        thisReader = FastAreader (self.infile)
        b1=int(self.bounds[0])
        b2=int(self.bounds[1])
        for head, seq in thisReader.readFasta():
            if b1<b2:
                print(f'>{head} ({self.bounds[0]},{self.bounds[1]}) (+)\n{seq[self.bounds[0]:self.bounds[1]+1]}\n')
            else:
                seq = self._reverseComplement(seq[self.bounds[1]:self.bounds[0]+1])
                print(f">{head} ({self.bounds[0]},{self.bounds[1]}) (5' --> 3') (-)\n{seq}\n")


    def _reverseComplement(self, seq):
        '''Returns the reverse complement in all-caps'''
        
        #found from www.educative.io, how to reverse a string
        complementDic = {'A':'T',
                         'T':'A',
                         'C':'G',
                         'G':'C'}
        rseq = ''.join(reversed(seq.upper()))
        rcseq = ''
        for nuc in rseq:
            rcseq += complementDic[nuc]
        return rcseq

def tester():
    loc = "C:/Users/mattk/Desktop/python scripts/scripting/Todd_intron_script/Todd_genomes/Acidovorax_avenae_subsp._avenae_ATCC_19860_wholegenome.fa"
    getBoundry = FastBoundry(bounds=(3980953,3980683), infile=loc)
    getBoundry.printBoundFasta()
#tester()
