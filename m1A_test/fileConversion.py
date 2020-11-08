from Bio import SeqIO

filepath = "TestingNNN/"
SeqIO.convert(filepath+'AlignedUsingEuk.sto', 'stockholm', filepath+'AlignedUsingEuk.fa', 'fasta')
print('hello')