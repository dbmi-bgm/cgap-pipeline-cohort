#!/bin/bash

SCRIPT_LOCATION="/usr/local/bin" # To use in prod

# variables from command line
input_vcf=$1
reference=$2
regionfile=$3
# data sources
vep_tar_gz=$4
dbnsfp_gz=$5
spliceai_snv_gz=$6
spliceai_indel_gz=$7
gnomad_gz=${8}
gnomad_gz2=${9}
CADD_snv=${10}
CADD_indel=${11}
# parameters
nthreads=${12}
version=${13} # 101
assembly=${14} # GRCh38

# self variables
directory=VCFS/

# rename with version
dbnsfp=dbNSFP4.1a.gz

# rename dbNSFP
ln -s $dbnsfp_gz $dbnsfp
ln -s ${dbnsfp_gz}.tbi ${dbnsfp}.tbi
ln -s ${dbnsfp_gz%.*}.readme.txt dbnsfp.readme.txt

# unpack data sources
tar -xzf $vep_tar_gz
tar -xzf ${vep_tar_gz%%.*}.plugins.tar.gz

# setting up output directory
mkdir -p $directory

# command line VEP
# plugins
plugin_dbnsfp="--plugin dbNSFP,${dbnsfp},GERP++_RS,GERP++_RS_rankscore,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,CADD_raw_rankscore,Polyphen2_HVAR_score,SIFT_pred,SIFT_converted_rankscore,SIFT_score,Ensembl_transcriptid"
plugin_spliceai="--plugin SpliceAI,snv=${spliceai_snv_gz},indel=${spliceai_indel_gz}"
plugin_CADD="--plugin CADD,${CADD_snv},${CADD_indel}"

plugins="--dir_plugins VEP_plugins --plugin TSSDistance $plugin_dbnsfp $plugin_spliceai $plugin_CADD"

# customs
custom_gnomad="--custom ${gnomad_gz},gnomADg,vcf,exact,0,AC,AF,AN"
custom_gnomad2="--custom ${gnomad_gz2},gnomADe2,vcf,exact,0,AC,AF,AN"

customs="$custom_gnomad $custom_gnomad2"

basic_vep="--sift b --polyphen b --symbol --canonical"

# options and full command line
options="--fasta $reference --assembly $assembly --use_given_ref --offline --cache_version $version --dir_cache . $basic_vep --force_overwrite --vcf --compress_output bgzip"


# Split file by chromosome and then in chunk of 500k variants
echo "Splitting files"
bcftools index -s $input_vcf | cut -f 1 > chromfile.txt
vep_chunk_file="./vep_chunk_files.txt"
command="bcftools view -O z --threads 8 -o split_by_chr.{}.vcf.gz $input_vcf {} || exit 1; python split_vcf.py -i $SCRIPT_LOCATION/split_by_chr.{}.vcf.gz -o $vep_chunk_file || exit 1; rm split_by_chr.{}.vcf.gz"
cat chromfile.txt | xargs -P $nthreads -i bash -c "$command" || exit 1 


# runnning VEP in parallel
echo "Running VEP"
command="vep -i {}.vcf.gz -o ${directory}{}.vep.vcf.gz $options $plugins $customs || exit 1; rm {}.vcf.gz || exit 1"
cat $vep_chunk_file | xargs -P $nthreads -i bash -c "$command" || exit 1

# merging the results
echo "Merging vcf.gz files"
array=(${directory}*.vep.vcf.gz)

IFS=$'\n' sorted=($(sort -V <<<"${array[*]}"))
unset IFS

files_sorted=""

for filename in ${sorted[@]};
  do
    #echo "Indexing file $filename"
    #tabix -p vcf -f "$filename" || exit 1
    files_sorted="$files_sorted$filename "
  done

echo "Concatenating files: $files_sorted"
bcftools concat -o combined.vep.vcf.gz --threads 24 -O z $files_sorted || exit 1
echo "Removing temporary files"
rm -f $files_sorted
# echo "Sorting and indexing combined file"
# bcftools sort -o combined.vep.vcf.gz -O z combined.vep.unsorted.vcf.gz || exit 1
echo "Indexing combined file"
tabix -p vcf combined.vep.vcf.gz || exit 1
