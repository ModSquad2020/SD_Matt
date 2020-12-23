import sys
sys.path.insert(1, "./Potentially massive folder")
sys.path.insert(2, "../../../CMPipelines/GenerateCM")
import os
from sequenceAnalysis import FastAreader
import pandas as pd

bigOlList = []

filepath = "./Potentially massive folder/speciesModDB/"
for org in os.listdir(filepath):
	#print(org)
	for mod in os.listdir(filepath+org):
		#print("\t"+mod)
		for pos in os.listdir(filepath+org+'/'+mod):
			#print("\t\t"+pos)
			fullFilePath = filepath+org+'/'+mod+'/'+pos
			for header, sequence in FastAreader(fullFilePath).readFasta():
					#print("\t\t\t"+header)
					genus, species, organismShort = org.split('_')[0], org.split('_')[1], org.split('_')[2]
					organismName = org.split('_')[0] + ' ' + org.split('_')[1]

					actualMod = mod.split('-')[len(mod.split('-'))-1]
					actualPos = pos.split('_')[ 1 ]
					
					isotype = header.split('-')[2]
					isoacceptor = header.split('-')[3]
					isodecoder = header.split('-')[4]
					bigOlList.append( [genus, species, organismName, organismShort, actualMod, actualPos, isotype, isoacceptor, isodecoder] )

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
df = pd.DataFrame(bigOlList, columns=['Genus', 'Species', 'FullName', 'ShortName', 'Modification', 'Position', 'Isotype', "Isoacceptor", "Isodecoder"])
df.to_csv('distilledModomics.csv', sep=' ', index=False)