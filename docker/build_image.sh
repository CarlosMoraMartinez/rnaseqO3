## Set up local registry
#sudo docker run -d -p 5000:5000 --restart always --name registry registry:2
# sudo vim /etc/docker/daemon.json
# { "insecure-registries": ["ccarlos.local:5000"] }

sudo docker build -t ccarlos/registry:samtools-1.20 -f samtools-1.20_image.dockerfile .
sudo docker build -t ccarlos/registry:STAR-2.7.11b -f star_image.dockerfile .
sudo docker build -t ccarlos/registry:HISAT2-2.2.1 -f hisat2-2.2.1_image.dockerfile .
sudo docker build -t ccarlos/registry:htslib-1.20 -f htslib-1.20_image.dockerfile .
sudo docker build -t ccarlos/registry:stringtie2-v2.2.2 -f stringtie2_image.dockerfile .
sudo docker build -t ccarlos/registry:salmondecoy -f salmontools_generate_decoy_image.dockerfile .
sudo docker build -t ccarlos/registry:kallisto-0.50.1 -f kallisto_image.dockerfile .

# Enter to test
sudo docker run -ti --entrypoint /bin/bash ccarlos/registry:kraken_with_pigz

# Enter and mount your home to the container
sudo docker run -ti  --mount type=bind,source=/home,target=/home  --entrypoint /bin/bash ccarlos/registry:HISAT2-2.2.1

# Save images in order to transfer to server
sudo docker save -o dockerimage_kraken_with_pigz.tar ccarlos/registry:kraken_with_pigz 
sudo docker save -o dockerimage_seqkit_with_samtools.tar ccarlos/registry:seqkit_with_samtools

# Transfer .tar files via scp and go to the server
cd /DATA12/COMUN/tmp
sudo docker load -i dockerimage_kraken_with_pigz.tar
sudo docker load -i dockerimage_seqkit_with_samtools.tar
