#!/bin/bash

# variables from command line
annotated_vcf=$1
annotated_vcf_tbi=$2
cases=$3
gene_annotations=$4
aaf_bin=$5
vc_tests=$6
excluded_genes=$7


# Optional arguments
if [ -z "$excluded_genes" ]
then
    excluded_genes="None"
fi

if [ -z "$vc_tests" ]
then
    vc_tests="burden"
fi


./run_regenie.sh -v $annotated_vcf \
                 -t $annotated_vcf_tbi \
                 -c $cases \
                 -g $gene_annotations \
                 -a $aaf_bin \
                 -b $vc_tests \
                 -e $excluded_genes || exit 1
