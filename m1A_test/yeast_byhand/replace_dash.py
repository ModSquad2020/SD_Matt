filename = '-trna_pos51_m1A_aligned.fa'
with open(filename, 'r') as file:
	with open('new_'+filename, 'w') as write_to:
		for line in file.readlines():
			line=line.replace('-','')
			write_to.write(line)
