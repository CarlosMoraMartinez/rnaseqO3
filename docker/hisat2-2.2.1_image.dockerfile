FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip && \
    apt-get install -y python3.10 && \
    ln -s /usr/bin/python3.10 /usr/bin/python3 && \
    ln -s /usr/bin/python3.10 /usr/bin/python



# From source
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download && \
    mkdir software && \
    mv download /software/hisat2-2.2.1.zip && \
    cd software && \
    unzip hisat2-2.2.1.zip && \
    rm hisat2-2.2.1.zip 

WORKDIR /software/hisat2-2.2.1/

RUN make
ENV PATH="${PATH}:/software/hisat2-2.2.1/"