
import sys
sys.path.insert(1, "C:/Users/mattk/Desktop/python scripts/todd stuff/Biotools/")
from sequenceAnalysis import FastAreader
from NNNSlicer import CreateWindows
import os
import copy

whereToMake='.\\'

'''
Commandline
python generateWindows.py fasta.fa secondarystructure.txt

* Will generate a directory called "fasta_Sliced"
* will populate the directory with fasta files
	named after their specified regions.
'''
fastafile = 'euk_sacc.fa'#sys.argv[1]
secondaryStructure = 'sacc_ss.txt'#sys.argv[2]


#important stuff
with open(secondaryStructure) as file:
	secondaryStructure = file.read()

def generateSub(original, collection, i):
	'''
	Given a combination and a deepcopied version, 
	wrapped in a for loop ranging from 1, 
	len(collection) for i
	'''
	base = True
	for item in collection:
		if len(item) > i:#pop the 2nd column
			base=False
			item.pop(i)

	cleanCollection = []
	for item in collection:#get rid of the top rows
		if len(item) <= i:
			pass
		else:
			cleanCollection.append(item)

	if base:
		return original
	else: 
		original.extend(copy.deepcopy( cleanCollection ))
		return generateSub( original, collection, i )
    
def cleanSub(original, culmination):
	for item in original:
		for item2 in culmination:
			if item == item2:
				culmination.remove(item)
	culmination.extend(original)

	return culmination

def generateCombinations():
    '''Returns a list of all combinations of the desiredRegions'''
    desiredRegions = list( CreateWindows().alignedRegions.keys() )
    
    combinations = []
    for i in range( len(desiredRegions)):
        division = []
        for j in range(i+1, len(desiredRegions)+1):
            goingForward = desiredRegions[i:j]
            division.append( goingForward )
        
        culmination=[]
        collection = division
        import copy
        for column in range(1, len(collection)):
            thing = copy.deepcopy(collection)
            culmination.extend(generateSub(copy.deepcopy(collection), thing, column) )
        combinations.extend( cleanSub(collection, culmination) )
        
    return combinations

def generateCombinations2():
	'''(DEPRECIATED) Returns a list of traditional combinations of the desiredRegions'''
	desiredRegions = list( CreateWindows().alignedRegions.keys() )
	print('start', desiredRegions)

	combinations = []
	for i in range( len(desiredRegions)+1 ):
		for j in range(i+1, len(desiredRegions)+1):
			goingForward = desiredRegions[i:j]
			goingBackward = desiredRegions[j:]
			goingBackward.insert(0, desiredRegions[i])
			print(goingBackward)

			if goingForward not in combinations:
				combinations.append( goingForward )
			if goingBackward and not (goingBackward in combinations):
				combinations.append( goingBackward )

	return combinations



def makeDirectory():
	directoryName = fastafile.split('.')[0]+'_Sliced'
	path = os.path.join(whereToMake, directoryName)
	try:
		os.mkdir(path)
	except FileExistsError:
		pass

	return path


#important stuff
#   0         1              2             3         4
#'AA-stem', 'D-arm', 'Anticodon loop', 'Introns', 'T-arm'
directory = makeDirectory()
for desiredRegions in generateCombinations():
	newFasta=''
	for header, sequence in FastAreader(fastafile).readFasta():
		windowMaker = CreateWindows(secondaryStructure, sequence=sequence, desiredRegions=desiredRegions)
		newFasta += f'>{header}\n{windowMaker.getWindowedSequence()}\n'
	
	fastaName = "_".join(desiredRegions) + '_Sliced.fa'
	with open( os.path.join(directory, fastaName), 'w' ) as file:
		file.write( newFasta )