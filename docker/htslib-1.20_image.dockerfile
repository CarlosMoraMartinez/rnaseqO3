FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip && \
    apt-get install -y python3.10 && \
    ln -s /usr/bin/python3.10 /usr/bin/python3 && \
    ln -s /usr/bin/python3.10 /usr/bin/python



# From source
RUN mkdir software && \
    cd software && \
    wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xvjf htslib-1.20.tar.bz2 && \
    rm htslib-1.20.tar.bz2


WORKDIR /software/htslib-1.20/
RUN ./configure --prefix=/usr/local/ && \
    make && \
    make install 

WORKDIR /

RUN rm -rf software