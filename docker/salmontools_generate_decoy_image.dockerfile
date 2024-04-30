FROM  ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip curl && \
    apt-get install -y wget make gcc && \
    apt-get install -y g++ && \
    apt-get install -y libz-dev xxd bzip2 libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev zlib1g libgsl-dev && \
    apt install -y software-properties-common lsb-release && \
    apt-get install -y build-essential libssl-dev && \
    apt-get install -y cmake && \
    apt-get install -y git 

# Bedtools From source
RUN mkdir software && \
    cd software && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xvf bedtools-2.31.1.tar.gz && \
    rm bedtools-2.31.1.tar.gz

WORKDIR /software/bedtools2/

RUN make

ENV PATH="${PATH}:/software/bedtools2/bin/"

# hstlib from source
WORKDIR /software/
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xvjf htslib-1.20.tar.bz2 && \
    rm htslib-1.20.tar.bz2


WORKDIR /software/htslib-1.20/
RUN ./configure --prefix=/usr/local/ && \
    make && \
    make install 

# Mashmap4 From source
WORKDIR /software/

RUN git clone https://github.com/marbl/MashMap

WORKDIR /software/MashMap

RUN cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build && \
    mv build/bin/* /usr/local/bin/

# Salmon tools
WORKDIR /software/

run git clone  https://github.com/COMBINE-lab/SalmonTools/ && \
    cp SalmonTools/scripts/* /usr/local/bin/  && \
    chmod +x /usr/local/bin/generateDecoyTranscriptome.sh  && \
    rm -rf 





