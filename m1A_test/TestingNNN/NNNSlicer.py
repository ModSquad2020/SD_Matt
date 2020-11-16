class CreateWindows:
    '''
    Given a secondary structure and a single sequence 
    this will convert it into just the regions we are interested in.
    
    Usage:
    for header, sequence in FastAreader(file.fa).readFasta():
        windowMaker = CreateWindows(secondaryStructure, sequence=sequence, desiredRegions=[])
        newFasta += f'>{header}\n{windowMaker.getWindowedSequence()}'
        
    Output:
    >Example
    GGGCGUGUGNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNN-----------NNNUCCGGUUCGAUUCCGGACUCGUCCA
    '''
    import re
    startAA = re.compile(r'\(+[^<>]{0,}')#r'\(+.{0,}?<'
    endAA = re.compile(r'\)+:')#r'.{1}\)+'
    loops = re.compile(r'<+.{3,}?>+')
    introns = re.compile(r'[^<>]*\.{3,}[^<>]*')
    patterns = [startAA, endAA, loops, introns]

    
    def __init__(self, secondaryStructure, sequence='', desiredRegions=['AA-stem','D-arm','Anticodon loop', 'T-arm', 'Introns']):
        self.ss = secondaryStructure
        self.alignedRegions = self.populateDictionary()
        if not self.isProperShape():
            raise Exception(self.ss, "Non-standard tRNA shape, or fix regex")
        self.userInput = sequence
        self.desiredRegions = desiredRegions
        
    def populateDictionary(self):
        patterns = self.patterns
        alignedRegions = {  'AA-stem': [],
                            'D-arm': [],
                            'Anticodon loop': [],
                            'Introns': [],
                            'T-arm': []
                         }
        
        for n, pattern in enumerate(patterns):
            #loops through the different regular expression patterns
            matches = pattern.finditer(self.ss)

            for matchNum, match in enumerate(matches):
                #loops through each match to the pattern
                if n==0 or n==1:
                    alignedRegions['AA-stem'].append( match.span() )
                elif n==2:
                    if matchNum==0:
                        alignedRegions['D-arm'].append( match.span() )
                    elif matchNum==1:
                        alignedRegions['Anticodon loop'].append( match.span() )
                    else:
                        alignedRegions['T-arm'].append( match.span() )
                else:
                    alignedRegions['Introns'].append( match.span() )
        
        return alignedRegions
    
    def isProperShape(self):
        '''Returns True if the dictionary has the proper shape of a tRNA'''
        ar = self.alignedRegions
        a = len( ar['AA-stem'] ) == 2
        b = len( ar['D-arm'] ) == 1
        c = len( ar['Anticodon loop'] ) == 1
        d = len( ar['T-arm'] ) == 1
        
        return a and b and c and d

    def converted2N(self):
        '''Returns a new fasta string, where every character is replaced by N'''
        newSequence = ''
        for n, letter in enumerate(self.userInput):
            if letter.isalpha():
                newSequence += "N"
            else:
                newSequence += "-"
        return newSequence
                
 
    def getDesiredRanges(self):
        '''Returns list of the ranges where we dont want just N'''
        ranges = []
        for rang in self.desiredRegions:
            ranges.extend( self.alignedRegions[rang] )
        return ranges
    
    def getWindowedSequence(self):
        newSequence = self.converted2N()
        for rang in self.getDesiredRanges():
            for position in range(rang[0], rang[1]):
                newSequence = newSequence[:position] + self.userInput[position] + newSequence[position+1:]
        return newSequence