from Bio import SeqIO

SeqIO.convert('sacCer3-trnaalign.stk', 'stockholm', 'sacCer3-trnaalign.fa', 'fasta')
print('hello')