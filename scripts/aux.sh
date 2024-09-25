

#Generate index STAR:

# Get fasta and gtf_
#  http://www.gencodegenes.org/mouse_releases/current.html
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir WhereYouWantIndex \
--genomeFastaFiles GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.vM11.primary_assembly.annotation.gtf --sjdbOverhang 100


# get fasta of transcripts
gffread -w mm39.knownGene.transcripts.fasta -g mm39.fa  genes/mm39.knownGene.gtf