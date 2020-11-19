from Bio import SeqIO
import sys

'''Converts stockholm --> fasta'''

print(sys.argv[1])
fromfile = sys.argv[1]
tofile = fromfile.split('.')[0]+'.fa'

SeqIO.convert(fromfile, 'stockholm', tofile, 'fasta')
print(f'Done: created {tofile}')
