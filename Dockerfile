FROM python:3

ADD bin/align_to_ref.py /usr/local/bin/align_to_ref.py
ADD bin/bam2consensus.py /usr/local/bin/bam2consensus.py

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN cp /minimap2-2.17_x64-linux/minimap2 /usr/local/bin

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2;tar -xvjf samtools-1.10.tar.bz2;cd samtools-1.10;./configure;make;make install

RUN pip install pysam biopython
RUN chmod +x /usr/local/bin/align_to_ref.py
RUN chmod +x /usr/local/bin/bam2consensus.py