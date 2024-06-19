FROM ubuntu:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y wget make gcc && \
    apt-get install -y g++ && \
    apt-get install -y libz-dev xxd bzip2 libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev && \
    apt-get install -y unzip && \
    apt-get install -y python3 && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    apt install -y python3-pip

RUN pip install deeptoolsintervals --break-system-packages && \
    pip install matplotlib --break-system-packages && \
    pip install numpydoc --break-system-packages && \
    pip install plotly --break-system-packages && \
    pip install py2bit --break-system-packages && \
    pip install pyBigWig --break-system-packages && \
    pip install scipy --break-system-packages && \
    pip install pysam --break-system-packages && \
    pip install HTSeq --break-system-packages 

