import matplotlib.pyplot as plt
import sys
sys.path.insert(1, 'C:\\Users\\mattk\\Desktop\\python scripts\\todd stuff\\BioTools')
from sequenceAnalysis import FastAreader
import math


NAMES = 'tRNAnames.txt'
EVALS = 'eValues.txt'

positiveFasta = "hasM1A.fa"

with open(NAMES) as nameFile:
	with open(EVALS ) as evalFile:
		NAMES = []
		EVALS = []
		COLOR = []
		for name, evalue in zip(nameFile.readlines(), evalFile.readlines()):
			#print(name[0:name.find(' ')], -math.log10(float( evalue )))
			NAMES.append(name[0:name.find(' ')])
			EVALS.append( -math.log10(float( evalue )) ) 

		headers=[]
		for header, sequence in FastAreader(positiveFasta).readFasta():
			headers.append(header.split(' ')[0])

		for name in NAMES:
			inTest = False
			for header in headers:
				#print(name, header)
				if name.find(header) > -1:
					inTest = True
					break

			if inTest:
				COLOR.append('green')
			else:
				COLOR.append('black')


		
		#print(len(COLOR), COLOR)
		x = NAMES
		x_pos = [i for i in x]

		y = EVALS

		plt.bar(x_pos, y, color = COLOR)
		plt.xticks(x_pos, x, rotation=90, fontsize=5)
		plt.tight_layout()

		plt.show()



