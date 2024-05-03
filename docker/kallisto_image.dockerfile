FROM ubuntu:latest

RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y wget make gcc && \
        apt-get install -y g++ && \
        apt-get install -y libz-dev xxd bzip2 libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev zlib1g libgsl-dev libhdf5-dev && \
        apt-get install -y cmake && \
        apt-get install -y git  

# Kallisto From source

RUN mkdir software

WORKDIR /software/

RUN git clone https://github.com/pachterlab/kallisto.git && \
    mkdir kallisto/build

WORKDIR /software/kallisto/build

RUN cmake .. && \
    make  && \
    make install

# Salmon tools
WORKDIR /software/

RUN rm -rf kallisto





