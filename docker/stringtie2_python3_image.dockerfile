FROM ccarlos/registry:samtools-1.20

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y unzip curl && \
    apt-get install -y git && \
    apt-get install -y python3 && \
    ln -s /usr/bin/python3 /usr/bin/python

# From source

RUN mkdir software

WORKDIR /software/

RUN git clone https://github.com/gpertea/stringtie

WORKDIR /software/stringtie/

RUN make release
ENV PATH="${PATH}:/software/stringtie/"

RUN wget https://raw.githubusercontent.com/CarlosMoraMartinez/rnaseqO3/feature/localcm/scripts/prepDE_mod.py3 && \
    chmod +x prepDE_mod.py3
