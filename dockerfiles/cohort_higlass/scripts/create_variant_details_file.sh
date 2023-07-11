#!/bin/bash
shopt -s extglob

echoerr() { 
    printf "%s\n" "$*" >&2
    exit 1
}

printHelpAndExit() {
    echo "Usage: ${0##*/} -a ANNOTATED_VCF -v HIGLASS_VCF -s SAMPLE_INFO"
    echo "-a VCF : path to VEP annotated VCF (gzipped)"
    echo "-s SAMPLE_INFO : JSON string with sample information"
    exit "$1"
}
while getopts "a:v:s:" opt; do
    case $opt in
        a) annotated_vcf="$OPTARG"
           annotated_tbi="$OPTARG.tbi"
        ;;
        s) sample_info=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "============================="
echo "Creating variant details file"
echo "============================="
echo "Annotated, filtered VCF: $annotated_vcf"
echo "Annotated index: $annotated_tbi"
echo ""
echo "Sample info: $sample_info" 
echo ""
echo "============================="


if [ -z "$annotated_vcf" ]
then
    echoerr "Annotated VCF missing"
fi

if [ -z "$annotated_tbi" ]
then
    echoerr "Annotated VCF index missing"
fi

if [ -z "$sample_info" ]
then
    echoerr "Sample info is missing"
fi


SCRIPT_LOCATION="/usr/local/bin" # To use in prod
#SCRIPT_LOCATION="/Users/alexandervelt/Documents/GitHub/cgap-pipeline-cohort/dockerfiles/regenie/scripts" # To use locally

echo ""
echo "== Create the file =="
python "$SCRIPT_LOCATION"/create_variant_details_file.py -a "$annotated_vcf" -s "$sample_info" -o variant_details.vcf.gz || exit 1

echo ""
echo "== DONE =="

