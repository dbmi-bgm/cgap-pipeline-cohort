FROM cgap/cgap-ubuntu2004-py-38:0.0.1
LABEL org.opencontainers.image.authors="Alexander Veit (alexander_veit@hms.harvard.edu)"

#######################################################################
#     Setting working env
#######################################################################
#/miniconda3/bin/regenie
WORKDIR /usr/local/bin

RUN apt-get update && apt-get install -y libzstd-dev \
    libboost-all-dev \
    locales \
    tabix \
    autoconf \
    libssl-dev \
    libbz2-dev \
    pkg-config \
    liblzma-dev

RUN ln -s /usr/bin/python3.8 /usr/bin/python


## vcftools
RUN git clone https://github.com/vcftools/vcftools.git && \
    mv vcftools vcftools_tmp && \
    cd vcftools_tmp && \
    git checkout 954e607 && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cd ..

## bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
    tar -xjf bcftools-1.11.tar.bz2 && \
    rm -rf bcftools-1.11.tar.bz2 && \
    cd bcftools-1.11 && \
    ./configure && \
    make && \
    cd .. && \
    ln -s bcftools-1.11 bcftools

RUN conda install -c bioconda  -c conda-forge \
    click \
    peddy=0.4.8


RUN pip install granite-suite==0.2.0 scipy

RUN conda clean -a -y -f



#######################################################################
#     Setting env variables
#######################################################################

ENV PATH=/usr/local/bin/bcftools/:$PATH
## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

#######################################################################
#     Scripts
#######################################################################


COPY scripts/apply_gatk_filter.py .
COPY scripts/filter_hwe_by_pop.pl .
RUN chmod +x filter_hwe_by_pop.pl
COPY scripts/create_hwe_popmap.py .
COPY scripts/run_peddy.py .
COPY scripts/run_filtering.sh .
RUN chmod +x run_filtering.sh
COPY scripts/split_multiallelics.sh .
RUN chmod +x split_multiallelics.sh


# default command
CMD ["run_filtering.sh"]
