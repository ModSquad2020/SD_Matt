import sys
sys.path.insert(1, "C:\\Users\\mattk\\Desktop\\python scripts\\todd stuff\\BioTools")
from sequenceAnalysis import FastAreader
import pandas as pd

'''
Given an aligned fasta (yeast) file
This program will write the output to a 
csv

'''

fileToDisplay = "sacCer3-trnaalign.fa"

tRNAdic = {}
for header, seq in FastAreader(fileToDisplay).readFasta():
	header = header.split(' ')[0]
	tRNAdic[header] = [position for position in seq]

#for yeast
columns = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
			'11', '12', '13', '14', '15', '16', '17a', '17b', '18', 
			'19', '20a', '20b', '20c', '21', '22', '23', '24', '25', 
			'26', '27', '28', '29', '30', '31', '32', '33', '34', 
			'35', '36', '37', '38', '39', '40', '41', '42', '43', 
			'44', '45', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 
			'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'e15', 
			'e16', 'e17', 'e18', 'e19', '46', '47', '48', '49', 
			'50', '51', '52', '53', '54', '55', '56', '57', '58', 
			'59', '60', '61', '62', '63', '64', '65a', '65b', '65c', 
			'65d', '65e', '70', '71', '72', '73', '74', '75', '76', 
           '77', '78', '79', '80']    
  
print(len(columns))
df = pd.DataFrame.from_dict(tRNAdic, orient='index', columns=columns)
pd.set_option('display.max_columns', None)
#pd.set_option('display.max_rows', None)
pd.options.display.max_colwidth = 1000
df.to_csv('sacCer_column.csv')
