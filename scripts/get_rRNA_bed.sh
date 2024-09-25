# conda install bioconda::bedops

gtf2bed < c_elegans.PRJNA13758.WS292.canonical_geneset.gtf | grep rRNA | grep -v MtDNA > c_elegans.PRJNA13758.WS292.rRNA.bed

bedtools merge -i c_elegans.PRJNA13758.WS292.rRNA.bed > c_elegans.PRJNA13758.WS292.rRNA.merged.bed

##

gtf2bed < c_elegans.PRJEB28388.WS292.annotations.gff3 | grep rRNA | grep -v MtDNA > c_elegans.PRJEB28388.WS292.rRNA.bed

bedtools merge -i c_elegans.PRJEB28388.WS292.rRNA.bed > c_elegans.PRJEB28388.WS292.rRNA.merged.bed

## Alternative way to get refflat (did not work sometimes):
gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons file.gtf tmp.txt
awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' tmp.txt > out.refflat
rm tmp.txt


# First download from BioMart, then get fasta of transcripts
gffread -w mm39.knownGene.transcripts.fasta -g mm39.fa  genes/mm39.knownGene.gtf