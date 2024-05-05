FROM ccarlos/registry:htslib-1.20

RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y wget make gcc && \
        apt-get install -y g++ && \
        apt-get install -y libz-dev xxd bzip2 libbz2-dev liblzma-dev zlib1g libgsl-dev libhdf5-dev autoconf && \
        apt-get install -y cmake && \
        apt-get install -y git  

# Kallisto From source

RUN mkdir software

WORKDIR /software/

RUN git clone https://github.com/pachterlab/kallisto.git && \
    mkdir kallisto/build

WORKDIR /software/kallisto/build

RUN cmake -DUSE_BAM=ON -DUSE_HDF5=ON .. && \
    make  && \
    make install

# Salmon tools
WORKDIR /software/

RUN rm -rf kallisto





