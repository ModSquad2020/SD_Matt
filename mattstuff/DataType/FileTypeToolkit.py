import pandas as pd
import numpy as np
from modomicsDictionary import modomicsDict as modDic

file = 'TESTresults.csv'
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

with open("modomics_organisms.txt") as file:
	with open("new_modomics_organisms.txt", 'w') as file2:
		newSet = set()
		for org in file.readlines():
			newSet.add(org)
		for item in newSet:
			file2.write(item)

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