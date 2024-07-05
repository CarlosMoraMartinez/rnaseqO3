
sudo docker pull staphb/fastqc:latest
#sudo docker pull bytesco/pigz:latest # Doesnt work
sudo docker pull nsheff/pigz:latest
sudo docker pull staphb/seqkit:latest
sudo docker pull staphb/multiqc:latest
sudo docker pull multiqc/multiqc
sudo docker pull pegi3s/multiqc:latest  #Use this one with version 1.14; the other ones with 1.22/23 do not recognise picard rnaseqmetrics

sudo docker pull staphb/trimmomatic:latest
sudo docker pull staphb/samtools:latest
sudo docker pull gcr.io/broad-cga-aarong-gtex/rnaseqc:latest

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

## own images at Dockerhub
docker pull carlosmora91/cmora_images:samtools-1.20
docker pull carlosmora91/cmora_images:STAR-2.7.11b
docker pull carlosmora91/cmora_images:htslib-1.20
docker pull carlosmora91/cmora_images:salmondecoy
docker pull carlosmora91/cmora_images:kallisto-0.50.1
docker pull carlosmora91/cmora_images:subread-2.0.6
docker pull carlosmora91/cmora_images:HISAT2-2.2.1
docker pull carlosmora91/cmora_images:stringtie2-v2.2.2
docker pull carlosmora91/cmora_images:HTSeq-2.0.5
docker pull carlosmora91/cmora_images:picard-3.1.1



