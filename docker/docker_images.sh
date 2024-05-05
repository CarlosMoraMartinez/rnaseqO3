
sudo docker pull staphb/fastqc:latest
#sudo docker pull bytesco/pigz:latest # Doesnt work
sudo docker pull nsheff/pigz:latest
sudo docker pull staphb/seqkit:latest
sudo docker pull staphb/multiqc:latest
sudo docker pull staphb/trimmomatic:latest
sudo docker pull staphb/samtools:latest

sudo docker pull pipecraft/cutadapt:0
sudo docker pull docker pull  
#sudo docker pull staphb/bbtools # Falla, parece que necesita instalar dependencias

sudo docker pull nanozoo/hisat2:2.1.0--66dae66 #2.1.0
sudo docker pull nanozoo/hisat2:2.2.1--75357cd # version 2.2.1

sudo docker pull mpgagebioinformatics/star:2.7.11b
# Install metaphlan database in permanent location:
# sudo docker run -ti --mount type=bind,source=/home,target=/home --entrypoint /bin/bash staphb/metaphlan
# metaphlan --install --bowtie2db /home/ccarlos/data/resources/metaphlan

sudo docker pull combinelab/salmon:latest