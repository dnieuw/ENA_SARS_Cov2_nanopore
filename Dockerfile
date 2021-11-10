FROM python:3.7

RUN apt-get update
RUN apt-get install -y tabix

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf - && \
    cp /minimap2-2.17_x64-linux/minimap2 /usr/local/bin
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar -xvjf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && make && make install

RUN pip install pysam biopython pymongo lxml

ADD bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*
