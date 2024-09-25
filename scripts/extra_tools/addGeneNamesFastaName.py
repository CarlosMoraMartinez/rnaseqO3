import os
import sys

# ARGS: fasta with cDNA

result = {}
ifasta = sys.argv[1]
ofasta = ifasta.replace('.fa', '_GeneNamesFirst.fa')
print(f"Output file: {ofasta}")

fo = open(ofasta, 'w')

for i in open(ifasta, 'r'):
    if i.startswith(">"):
        if 'gene_symbol:' in i:
            i = f">{[j.split(':')[1] for j in i.split(' ') if j.startswith('gene_symbol')].pop()}.{i[1::]}"
    fo.write(i + '\n')
       
fo.close()

   
