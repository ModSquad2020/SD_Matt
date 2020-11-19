import sys
sys.path.insert(1, "C:/Users/mattk/Desktop/python scripts/todd stuff/Biotools/")
from sequenceAnalysis import FastAreader
import os

'''
Commandline usage:
python fastaFileSorter.py "positiveModomics.fa" "unsorted.fa"

output:
2 files in the same directory as this program
	* unsorted_pos.fa
	* unsorted_neg.fa

Notes:
The headers in the positiveModomics.fa have to be in the Henry format
The headers in teh unsorted.fa have to be in the gtRNAdb 11/18/2020 format
'''


sortedModomics = sys.argv[1]
modifier = sortedModomics.replace('-','_').split('_')[2]

fastaFile = sys.argv[2]
base,ext = fastaFile.split('.')

newFileName_pos = f'{base}_{modifier}.{ext}'
newFileName_neg = f'{base}_neg.{ext}'

sortedHeaderSet = set()
headerSet = set()
for sortedHeader, sortedSequence in FastAreader(sortedModomics).readFasta():
	sortedHeaderSet.add( tuple(sortedHeader.replace('Ini', 'iMet').split('-')[2:4]) )
	pass

with open(newFileName_pos, 'w') as posFile:
	with open(newFileName_neg,'w') as negFile:
		for header, sequence in FastAreader(fastaFile).readFasta():
			thing = tuple(header.split('-')[1:3])
			headerSet.add( tuple(header.split('-')[1:3]) )
			if thing in sortedHeaderSet:
				posFile.write( f'>{header}\n{sequence}\n' )
			else:
				negFile.write( f'>{header}\n{sequence}\n' )


	

