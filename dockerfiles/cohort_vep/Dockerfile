#######################################################################
#     Basic image
#######################################################################
FROM cgap/cgap-ubuntu2004-py-38:0.0.1
LABEL org.opencontainers.image.authors="Alexander Veit (alexander_veit@hms.harvard.edu)"

#######################################################################
#     General updates & installing necessary Linux components
#     VEP:
#       - libbz2-dev, liblzma-dev (HTSlib)
#       - libdbi-perl, libmodule-build-perl,
#         libarchive-zip-perl, libtest-exception-perl (Perl)
#       - libmysqlclient-dev, libpng-dev,
#         libexpat1-dev, libssl-dev (other)
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix \
    libbz2-dev \
    liblzma-dev \
    libdbi-perl \
    libmodule-build-perl \
    libarchive-zip-perl \
    libtest-exception-perl \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    libxml-libxml-perl \
    libexpat1-dev \
    autoconf \
    pkg-config 

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

#######################################################################
#     Software
#######################################################################


## VEP
RUN git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    git checkout release/101 && \
    perl INSTALL.pl --AUTO a --NO_UPDATE && \
    cd ..

## perl module for VEP
RUN cpan App::cpanminus && \
    wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz && \
    tar xzf v335_base.tar.gz && \
    export KENT_SRC=/usr/local/bin/kent-335_base/src && \
    export MACHTYPE=$(uname -m) && \
    export CFLAGS="-fPIC" && \
    MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'` && \
    export MYSQLINC && \
    MYSQLLIBS=`mysql_config --libs` && \
    export MYSQLLIBS && \
    cd $KENT_SRC/lib && \
    echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk && \
    make clean && make && \
    cd ../jkOwnLib && \
    make clean && make && \
    cd /usr/local/bin/ && \
    mkdir -p /usr/local/bin/cpanm_bw && \
    export PERL5LIB=$PERL5LIB:/usr/local/bin/cpanm_bw/lib/perl5 && \
    cpanm -l /usr/local/bin/cpanm_bw Bio::DB::BigFile

## Test
RUN export PERL5LIB=$PERL5LIB:/usr/local/bin/ensembl-vep && \
    cd /usr/local/bin/ensembl-vep && \
    perl -Imodules t/AnnotationSource_File_BigWig.t && \
    cd /usr/local/bin


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

RUN pip install granite-suite==0.2.0 click

#######################################################################
#     Scripts
#######################################################################
## vep-annot
COPY vep-annot.sh .
RUN chmod +x vep-annot.sh
COPY split_vcf.py .

#######################################################################
#     Setting env variables
#######################################################################
ENV PATH=/usr/local/bin/bcftools/:$PATH
ENV PATH=/usr/local/bin/ensembl-vep/:$PATH
ENV PERL5LIB=$PERL5LIB:/usr/local/bin/cpanm_bw/lib/perl5:/usr/local/bin/ensembl-vep

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
