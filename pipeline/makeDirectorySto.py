'''
Usage:
python makeDirectorySto.py "directoryName" "secondaryStructure.txt"

Output:
A directory of .sto files

Notes:
	* This directory needs to only be full of .fa files
	* This stockholm is different from normal, but it should build the same
'''

'''
Given a directory name and a text file with the secondary structure
convert all the .fa files into .sto files
'''
from Bio import SeqIO
import sys
import os



directoryName = sys.argv[1]
ssFile = sys.argv[2]

ss = ''
with open(ssFile) as file:
	ss = file.read()

for filename in os.listdir(directoryName):
	path = os.path.join(directoryName, filename)

	tofile = filename.split('.')[0]+'.sto'
	path2 = os.path.join(directoryName, tofile)

	SeqIO.convert(path, 'fasta', path2, 'stockholm')
	
	print(filename, '   ', tofile)

	newContents = ''
	with open(path2, 'r') as file:
		contents = file.read()
		toReplace = f'#=GC SS_cons {ss}//'

		newContents = contents.replace('//', toReplace)
	with open(path2, 'w') as file:
		file.write( newContents )
	
	os.remove(path)
	