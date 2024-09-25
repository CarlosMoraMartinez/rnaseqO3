
result = []
pattern = "Il22"
getseq = False

for i in open('Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all.fa', 'r'):
    if i.startswith(">"):
        if pattern in i:
            getseq = True
            result.append(i.strip())
        else:
            getseq = False
    elif getseq:
        result.append(i.strip())
       
fo = open(f"{pattern}_filt.fa", 'w')

for l in result:
    fo.write(l + '\n')
    
fo.close()

   
