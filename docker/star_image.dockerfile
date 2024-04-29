FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
        apt-get upgrade -y

# Install STAR aligner
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz  && \
    tar -xzf 2.7.11b.tar.gz 

WORKDIR STAR-2.7.11b/source

RUN    make STAR && \
    mv STAR /usr/local/bin

WORKDIR /

RUN rm -rf STAR-2.7.11b STAR-2.7.11b.tar.gz