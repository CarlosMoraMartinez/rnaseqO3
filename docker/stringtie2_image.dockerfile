FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip curl && \
    apt-get install -y git && \
    apt-get install -y python2.7 && \
    ln -s /usr/bin/python2.7 /usr/bin/python2 && \
    ln -s /usr/bin/python2.7 /usr/bin/python

# Note: the prepDE.py script to merge stringtie results uses python 2

# From source
RUN mkdir software && \
    cd software && \
    git clone https://github.com/gpertea/stringtie


WORKDIR /software/stringtie/

RUN make release
ENV PATH="${PATH}:/software/stringtie/"