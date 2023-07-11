#!/bin/bash
shopt -s extglob

echoerr() { 
    printf "%s\n" "$*" >&2
    exit 1
}

printHelpAndExit() {
    echo "Usage: ${0##*/} -v VCF -s SAMPLE_INFO -g GENE_ANNOTATIONS -a AAF_BIN"
    echo "-v VCF : path to VEP annotated VCF (gzipped)"
    echo "-s SAMPLE_INFO : JSON string with sample information"
    echo "-a AAF_BIN : specifies the AAF upper bound used to generate burden masks"
    echo "-r AF_THRESHOLD_HIGLASS : AF threshold for rare variants that are included in Higlass (e.g. 0.03)"
    echo "-g GENE_ANNOTATIONS : gene annotation file from portal"
    exit "$1"
}
while getopts "v:s:g:a:b:c:r:e:" opt; do
    case $opt in
        v) annotated_vcf="$OPTARG"
           annotated_vcf_tbi="$OPTARG.tbi"
        ;;
        s) sample_info=$OPTARG;;
        g) gene_annotations=$OPTARG;;
        a) aaf_bin=$OPTARG;;
        r) af_threshold_higlass=$OPTARG;;
        b) regenie_variant_results=$OPTARG;;
        c) regenie_gene_results=$OPTARG;;
        d) regenie_gene_results_snplist=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "============================="
echo "Annotated VCF: $annotated_vcf"
echo "Annotated VCF index: $annotated_vcf_tbi"
echo "Gene annotation file: $gene_annotations"
echo ""
echo "Sample info: $sample_info"
echo ""
echo "============================="

if [ -z "$annotated_vcf" ]
then
    echoerr "Annotated VCF missing"
fi

if [ -z "$annotated_vcf_tbi" ]
then
    echoerr "Annotated VCF index missing"
fi

if [ -z "$sample_info" ]
then
    echoerr "Sample info is missing"
fi

if [ -z "$gene_annotations" ]
then
    echoerr "Gene annotations missing"
fi

if [ -z "$aaf_bin" ]
then
    echoerr "AAF_BIN missing"
fi

if [ -z "$af_threshold_higlass" ]
then
    echoerr "Rare variant AF threshold for Higass result files missing."
fi

if [ -z "$regenie_variant_results" ]
then
    echoerr "Regenie variants results missing"
fi

if [ -z "$regenie_gene_results" ]
then
    echoerr "Regenie gene results missing"
fi

if [ -z "$regenie_gene_results_snplist" ]
then
    echoerr "Regenie gene snplist missing"
fi


SCRIPT_LOCATION="/usr/local/bin" # To use in prod
#SCRIPT_LOCATION="/Users/alexandervelt/Documents/GitHub/cgap-pipeline-cohort/dockerfiles/regenie/scripts" # To use locally


# STOPPING HERE FOR TESTING
echo 'This is a test' > variant_level_results.txt.gz
echo 'This is a test' > higlass_variant_tests.multires.vcf.gz
echo 'This is a test' > higlass_variant_tests.multires.vcf.gz.tbi
echo 'This is a test' > higlass_gene_tests.sorted.vcf.gz
echo 'This is a test' > higlass_gene_tests.sorted.vcf.gz.tbi
echo 'This is a test' > coverage.bw
exit 0

echo ""
echo "== Create coverage bigWig file =="
create-coverage-bed -i "$annotated_vcf" \
                  -o coverage.bed \
                  -a hg38 \
                  -q False || exit 1

convert-bed-to-bw -i coverage.bed \
                  -o coverage.bw \
                  -a hg38 \
                  -l 0 || exit 1



echo ""
echo "== Create variant level result file and Higlass VCF =="

python "$SCRIPT_LOCATION"/create_variant_result_file.py -r "$regenie_variant_results" \
                                      -a "$annotated_vcf" \
                                      -s "$sample_info" \
                                      -o variant_level_results.txt.gz \
                                      -f "$af_threshold_higlass" \
                                      -e higlass_variant_tests.gz || exit 1

# higlass_variant_tests.gz is gzip compressed. Recompress here with bgzip
gzip -cd higlass_variant_tests.gz | bgzip --threads 6 -c > higlass_variant_tests.vcf.gz || exit 1
rm -f higlass_variant_tests.gz
tabix -p vcf higlass_variant_tests.vcf.gz || exit 1

echo ""
echo "== Create multilevel version of the Higlass VCF =="
# This needs "pip install cgap-higlass-data"
# Skip sorting for now. It should already be sorted
#cat higlass_variant_tests.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_variant_tests.sorted.vcf

# Output will be compressed and indexed
create-cohort-vcf -i higlass_variant_tests.vcf.gz \
                  -o higlass_variant_tests.multires.vcf.gz \
                  -c fisher_ml10p_control \
                  -q True \
                  -t True \ 
                  -w True || exit 1

# cat higlass_variant_tests.multires.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_variant_tests.multires.sorted.vcf || exit 1
# rm -f higlass_variant_tests.multires.vcf
# bgzip -c higlass_variant_tests.multires.sorted.vcf > higlass_variant_tests.multires.vcf.gz || exit 1
# rm -f higlass_variant_tests.multires.sorted.vcf
# tabix -p vcf higlass_variant_tests.multires.vcf.gz || exit 1


echo ""
echo "== Create gene level Higlass file =="

# Extract gene annotations if necessary
if file --mime-type "$gene_annotations" | grep -q gzip$; then
  gunzip -c "$gene_annotations" > gene_annotations.tsv
else
  cp "$gene_annotations" gene_annotations.tsv
fi

python "$SCRIPT_LOCATION"/create_higlass_gene_file.py -r "$regenie_gene_results" \
                                   -g gene_annotations.tsv \
                                   -s "$regenie_gene_results_snplist" \
                                   -a "$aaf_bin" \
                                   -o higlass_gene_tests.vcf || exit 1

cat higlass_gene_tests.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_gene_tests.sorted.vcf || exit 1
bgzip -c higlass_gene_tests.sorted.vcf > higlass_gene_tests.sorted.vcf.gz || exit 1
tabix -p vcf higlass_gene_tests.sorted.vcf.gz || exit 1


echo ""
echo "== DONE =="

