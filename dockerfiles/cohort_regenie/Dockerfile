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
    gzip \
    less

RUN ln -s /usr/bin/python3.8 /usr/bin/python

RUN conda install -c bioconda  -c conda-forge \
    regenie=3.2.7 \
    plink2 \
    vcftools \
    bgenix \
    click 


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

COPY scripts/utils.py .
COPY scripts/create_mask_files.py .
COPY scripts/create_phenotype.py .
COPY scripts/run_regenie.sh .
RUN chmod +x run_regenie.sh

# default command
CMD ["run_regenie.sh"]
