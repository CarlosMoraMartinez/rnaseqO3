import os
import sys

cdna_fasta = sys.argv[1]
genelist_file = sys.argv[2]
outdir = sys.argv[3]

genelist = [i.strip() for i in open(genelist_file, 'r')]


result = {}
#pattern = "Il22"
getseq = False

for i in open(cdna_fasta, 'r'):
    if i.startswith(">"):
        for pattern in genelist:
            if f"gene_symbol:{pattern} " in i or f"gene:{pattern}" in i:
                print(f"gene {pattern} found in sequence {i}")
                getseq = True
                result.setdefault(pattern, [])
                result[pattern].append(i.strip())
                break
        else:
            getseq = False
    elif getseq:
        result[pattern].append(i.strip())
       
try:
    os.system(f"mkdir {outdir}")
except:
    print("Output dir already exists")

for k, v in result.items():
    fo = open(f"{outdir}/{k}_filt.fa", 'w')
    for l in v:
        fo.write(l + '\n')    
    fo.close()

   
