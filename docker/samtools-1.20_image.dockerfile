FROM ubuntu:latest

RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y wget make gcc && \
        apt-get install -y g++ && \
        apt-get install -y libz-dev xxd bzip2 libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev


# Instlal samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvjf samtools-1.20.tar.bz2

WORKDIR samtools-1.20

RUN ./configure --prefix=/usr/local/ && \
    make && \
    make install 

WORKDIR /

RUN rm -rf samtools-1.20*