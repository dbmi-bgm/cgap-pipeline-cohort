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
    less

RUN ln -s /usr/bin/python3.8 /usr/bin/python

RUN conda install -c bioconda  -c conda-forge \
    vcftools \
    bcftools>=1.11 \
    click \
    peddy=0.4.8



RUN pip install granite-suite==0.2.0 scipy

RUN conda clean -a -y -f



#######################################################################
#     Setting env variables
#######################################################################

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
COPY split_multiallelics.sh .
RUN chmod +x split_multiallelics.sh


# default command
CMD ["run_filtering.sh"]
