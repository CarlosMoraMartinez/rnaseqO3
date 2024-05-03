FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip curl

RUN mkdir software 

WORKDIR /software

# subread from source

RUN wget https://downloads.sourceforge.net/project/subread/subread-2.0.6/subread-2.0.6-source.tar.gz && \
    tar -xvf subread-2.0.6-source.tar.gz && \
    rm subread-2.0.6-source.tar.gz

WORKDIR /software/subread-2.0.6-source/src

RUN make -f Makefile.Linux

ENV PATH="${PATH}:/software/subread-2.0.6-source/bin/"

WORKDIR /