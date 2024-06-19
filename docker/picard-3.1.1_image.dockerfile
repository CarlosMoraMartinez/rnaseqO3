FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y libcurl4-openssl-dev


# Install JAVA

RUN mkdir software

WORKDIR /software

RUN wget https://download.oracle.com/java/22/latest/jdk-22_linux-x64_bin.deb && \
    mkdir /usr/share/binfmts/ && \
    dpkg -i jdk-22_linux-x64_bin.deb  && \
    rm jdk-22_linux-x64_bin.deb

# Install picard

RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar

ENV PICARD="/software/picard.jar"

# Install gtfToGenePred to create refflat
# For all UCSC tools:
# rsync -aP   hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./

WORKDIR /software

RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && \
    mv gtfToGenePred /usr/local/bin && \
    chmod a+x /usr/local/bin/gtfToGenePred

# Download GATK

RUN apt-get install unzip && \
    wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip gatk-4.5.0.0.zip && \
    rm gatk-4.5.0.0.zip 

ENV GATK="/software/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"

# Download GTF2refflat

RUN  wget https://github.com/biopet/gtftorefflat/releases/download/v0.1/GtftoRefflat-assembly-0.1.jar

ENV GTFTOREFFLAT="/software/GtftoRefflat-assembly-0.1.jar"

WORKDIR /