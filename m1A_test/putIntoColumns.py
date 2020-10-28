import sys
sys.path.insert(1, "C:\\Users\\mattk\\Desktop\\python scripts\\todd stuff\\BioTools")
from sequenceAnalysis import FastAreader
import pandas as pd

tRNAdic = {}
for header, seq in FastAreader('alignedYeasttRNA.fa').readFasta():
	header = header.split(' ')[0]
	tRNAdic[header] = [position for position in seq]

df = pd.DataFrame.from_dict(tRNAdic, orient='index')
pd.set_option('display.max_columns', None)
print(df)
