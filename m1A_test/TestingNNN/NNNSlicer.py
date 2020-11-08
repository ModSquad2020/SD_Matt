'''


|       AA-acceptor          ||                D-arm                       ||                                          Anitcodon Loop                                                         ||                        Variable Loop              ||                    T-arm                         ||       AA-Acceptor         |
(.(.(..(....(..(.(......,..,..<.<<.<._......__...._....___.........._>.>.>.>....,.<....<..<.<.<____.__............................................................................._.>.>..>.>.>,...,<<..<<<<<.....___.............._..>>>>>.>..>,..,<.<...<..<.<._...._...._...._._......_.._..>.>>>...>......).).)...).)...)...)...:
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234
|0       |1        |2        |3        |4        |5        |6        |7        |8        |9        |10       |11       |12       |13       |14       |15       |16       |17       |18       |19       |20       |21       |22       |23       |24       |25       |26       |27       |28       |29       |30       |31       |32             

=============================================================

|   aa1 ||    D-arm       ||    AC    ||                        Intron                            || AC || Variable loop ||    T-arm      ||  aa2 |
(((((((,,<<<<___.___.__>>>>,<<<<<______............................................................_>>>>>,,............,,,<<<<<_______>>>>>))))))):
012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
|0       |1        |2        |3        |4        |5        |6        |7        |8        |9        |10       |11       |12       |13       |14    

'''

#alignedagainstEU.sto secondary strucutre
alignedRegions = {  'AA-acceptor1': (0,9),
					'AA-acceptor2': (139,147),
					'D-arm': (9,27),
					'Anticodon loop1': (27,39),
					'Intron': (39,99),
					'Anticodon loop2': (99,105),
					'Variable loop': (105,122),
					'T-arm': (122,139)

				}

def oneRegion(primarySequence, na):
	part1=na1
	a1= alignedRegions[part1][0]
	a2= alignedRegions[part1][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(a1,a2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence

def twoRegions(primarySequence, na1, na2):
	part1=na1
	a1= alignedRegions[part1][0]
	a2= alignedRegions[part1][1]

	part2=na2
	b1= alignedRegions[part2][0]
	b2= alignedRegions[part2][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(a1,a2) or n in range(b1,b2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence

def threeRegions(primarySequence, na1, na2, na3):
	''''''
	part1=na1
	a1= alignedRegions[part1][0]
	a2= alignedRegions[part1][1]

	part2=na2
	b1= alignedRegions[part2][0]
	b2= alignedRegions[part2][1]

	part3=na3
	c1= alignedRegions[part3][0]
	c2= alignedRegions[part3][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(a1,a2) or n in range(b1,b2) or n in range(c1,c2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence

'''
#if i am not lazy implement these
def fourRegions(na1, na2, na3, na4):
	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(f1,f2) or n in range(s1,s2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence
def fiveRegions(na1, na2, na3, na4, na5):
	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(f1,f2) or n in range(s1,s2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence
def sixRegions(na1, na2, na3, na4, na5, na6):
	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(f1,f2) or n in range(s1,s2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence
def sevenRegions(na1, na2, na3, na4, na5, na6, na7):
	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	part1=na
	a1= alignedRegions[part2][0]
	a2= alignedRegions[part2][1]

	newSequence = ""
	for n, letter in enumerate(primarySequence):
		if letter.isalpha():
			if n in range(f1,f2) or n in range(s1,s2):
				newSequence += letter
			else:
				newSequence += "N"
		else:
			newSequence += "-"
	return newSequence

'''

'''
#yeaststuff.sto secondary structure
alignedRegions = {  'AA-acceptor1': (0,30),
					'AA-acceptor2': (296,325),
					'D-arm': (30,76),
					'Anticodon loop': (76,191),
					'Variable loop': (191,244),
					'T-arm': (244,296)

				}

secondaryStructure="(.(.(..(....(..(.(......,..,..<.<<.<._......__...._....___.........._\
>.>.>.>....,.<....<..<.<.<____.__.....................................\
........................................_.>.>..>.>.>,...,<<..<<<<<....\
.___.............._..>>>>>.>..>,..,<.<...<..<.<._...._...._...._._.....\
._.._..>.>>>...>......).).)...).)...)...)...:"

part='D-arm'
darm=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]
part='T-arm'
tarm=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]
part='Variable loop'
vloop=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]
part='Anticodon loop'
antiloop=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]
part='AA-acceptor1'
acceptor1=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]
part='AA-acceptor2'
acceptor2=secondaryStructure[alignedRegions[part][0]: alignedRegions[part][1]]

parts=[darm,tarm,vloop,antiloop,acceptor1,acceptor2]


for part in parts:
	print(len(part), part)

print( sum(len(part) for part in parts) )
'''