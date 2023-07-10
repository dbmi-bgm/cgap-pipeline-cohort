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
    vcftools \
    click \
    ucsc-bedgraphtobigwig

# RUN conda install -c bioconda  -c conda-forge \
#     peddy=0.4.8 \
#     ucsc-bedgraphtobigwig

RUN pip install cgap-higlass-data==0.3.0 granite-suite==0.2.0 scipy

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
COPY scripts/create_higlass_gene_file.py .
COPY scripts/create_variant_result_file.py .
# COPY scripts/run_peddy.py .
COPY scripts/gather_results.sh .
RUN chmod +x gather_results.sh
COPY scripts/create_variant_details_file.py .
COPY scripts/create_variant_details_file.sh .
RUN chmod +x create_variant_details_file.sh

# default command
CMD ["gather_results.sh"]
