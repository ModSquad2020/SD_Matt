import pandas as pd
import numpy as np
import sys
sys.path.insert(1, "./Potentially massive folder")
from modomicsDictionary import modomicsDict as modDic
import csv
from copy import deepcopy

file = 'TESTresults.csv'
file2 = "distilledModomics.csv"
delim = ' '

def addEmptyColumn(name):
	df = pd.read_csv(file, delimiter=delim)
	if name not in list(df.columns):
		df[name] = np.nan
		df.to_csv(file, sep=delim, index=False)

def removeColumn(name):
	df = pd.read_csv(file, delimiter=delim)
	if name in list(df.columns):
		df.drop(columns=name, inplace=True)
		df.to_csv(file, sep=delim, index=False)

#listOfMods = ''
'''
with open('mod_shortnames.txt', 'w') as file:
	for key, values in modDic.items():
		print(values['short_name'])
		thing=values['short_name']
		if thing not in 'AGCU':
			try:
				newThing = thing#.replace('(', "\\(").replace(')', "\\)")
				listOfMods += f' {newThing}'
				file.write(f'{thing}\n')
			except UnicodeEncodeError:
				pass
	print(listOfMods)

'''
#for translating the dumb modomics format anticodons into cannonical bases
originalBaseDic = {'τ': 'C', 'F': 'U', '∨': 'G', '∞': 'A', 'ς': 'G', '4': 'U', 
					'α': 'U', 'ν': 'C', '`': 'A', '℘': 'U', 'T': 'U', 'R': 'G', 'D': 'U', 'μ': 'C', 
					'E': 'A', 'ρ': 'U', '?': 'C', '√': 'A', '^': 'A', 'Ê': 'U', 'K': 'G', ')': 'U', 
					'Î': '', 'O': 'A', '1': 'U', '#': 'G', '~': 'U', 'Y': 'G', '>': 'C', 'œ': 'A', 
					'ξ': 'A', 'ℑ': 'G', '∫': 'U', 'L': 'G', '/': 'A', 'G': 'G', 'λ': 'C', '|': 'G', 
					'ζ': 'A', '⊄': 'G', '&': 'U', '◊': 'U', '$': 'U', 'C': 'C', 'V': 'U', 'κ': 'U', 
					'δ': 'U', 'J': 'U', '*': 'A', "'": 'C', '⊥': 'U', '[': 'A', '+': 'A', 'σ': 'U', 
					'υ': 'U', '∑': 'G', 'U': 'U', '∇': 'G', 'π': 'U', 'ε': 'G', '∝': 'U', 'Q': 'G', 
					';': 'G', '3': 'U', 'β': 'C', '<': 'C', '¥': 'G', 'H': 'A', '≥': 'U', '5': 'U', 
					'X': 'U', 'N': 'U', 'Ω': 'G', '≅': 'U', '≈': 'A', '∩': 'U', 'M': 'C', 'χ': 'A', 
					'"': 'A', 'B': 'C', ',': 'U', '9': 'G', '!': 'U', ']': 'U', '⊇': 'G', '8': 'G', 
					'∃': 'U', '%': 'C', 'š': 'G', 'ω': 'U', 'ψ': 'G', '=': 'A', '⊆': 'G', '@': '', 
					'(': 'G', '2': 'U', '}': 'C', 'γ': 'G', 'Z': 'U', '°': 'C', '∴': 'G', 'ℵ': 'C', 
					'†': 'G', 'A': 'A', 'S': 'U', '∠': 'G', '⇑': 'G', 'I': 'A', '≡': 'A', '∪': 'U', 
					'7': 'G', '∏': 'U', 'P': 'U', '∅': 'C', 'φ': 'G', 'W': 'G', '\\': 'U', '∉': 'G', 
					'η': 'A', '⇓': 'A', '≠': 'A', ':': 'A', '≤': 'A', '{': 'U', '6': 'A', 'â': 'A', 
					'©': 'pppN', '®': 'pppN', '¶': 'pppN', '§': 'pppN', '±': 'A', 'Ð': 'U', 'Þ': 'U', 
					'€': 'G', '.': '', '¿': 'C', 'æ': 'G', 'Ç': 'C', '¾': 'U', '¼': 'U', '½': 'U', 
					'ϖ': 'pppN', 'ϑ': 'pppN', '♠': 'pppN', '♣': 'pppN', '♥': 'pppN', '♦': 'pppN', 
					'ϒ': 'pppN', 'Ξ': 'pppN', 'y': 'G', 'none': 'pppN', 'l': 'U', 'r': 'U', 'b': 'U', 
					'm': 'pppN', 'Γ': 'U', 'f': 'U', 'h': 'U', 'Δ': 'U', '': 'pppN', 'e': 'A', 'Ϫ': 'A', 
					'Ϩ': 'A', 'Ѷ': 'U', 'Ͽ': 'U', 'ÿ': 'A', '«': 'A', '£': 'A', '¡': 'C'}


andrewList = []
with open(file) as csvfile:
	csvfile.readline()
	andrewReader = csv.reader(csvfile, delimiter = delim)
	for line in andrewReader:
		andrewList.append(line)

distilledModomics = []
with open(file2) as csvfile:
	csvfile.readline()
	andrewReader = csv.reader(csvfile, delimiter = delim)
	for line in andrewReader:
		distilledModomics.append(line)


num = 0
newThing = deepcopy(andrewList)
print(distilledModomics[0][9])

for andrewItem in andrewList:
	for modItem in distilledModomics:
		if modItem[6] == 'Ini':
			modItem[6] = 'iMet'

		if andrewItem[0] == modItem[2]:
			#if names are the same
			if f'{modItem[6]}-{modItem[7]}-{modItem[8]}' == '-'.join( andrewItem[1].split('-')[0:3] ):
				#if tRNAs are the same
				if andrewItem[2].split('_')[2] == modItem[4]:
					#if mod is the same
					if andrewItem[2].split('_')[3] == modItem[5]:
						#if position is the same
						toAdd = [modItem[9], '+']
						andrewItem.extend(toAdd)
						break
	try:
		test = andrewItem[10] + ""
	except IndexError:
		toAdd = [None, '-']
		andrewItem.extend(toAdd)
'''
for n, thing in enumerate( andrewList ):
	print(n, thing, '\n')
'''

df = pd.DataFrame( andrewList, columns = ['Organism Name', 'tRNA', "Query CM", 'GC Content', "Description", "pos Score", "pos eVal", "neg Score", "neg eVal", "Modomics Sequence", "Modomcis Positive"] )
df.to_csv("isItTheSame.csv", sep=delim, index=False)
		

	