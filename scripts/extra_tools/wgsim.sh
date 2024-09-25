
FASTA='../Mus_musculus_c57bl6nj.C57BL_6NJ_v1/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all_GeneNamesFirst.fa'
OUTDIR=simulated_fastq_bygene_error0
ONAME=$OUTDIR'/cDNAall_'

wgsim -e 0.0 -d 300 -s 50 -N 50000000 -1 150 -2 150 -r 0.001 -R 0 -X 0 -S 123 $FASTA $ONAME'R1.fastq' $ONAME'R2.fastq' >$ONAME'WGSIM.log' 2>$ONAME'WGSIM.err'; 

pigz $ONAME'R1.fastq'
pigz $ONAME'R2.fastq'