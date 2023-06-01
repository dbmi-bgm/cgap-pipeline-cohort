#!/bin/bash
shopt -s extglob

echoerr() { 
    printf "%s\n" "$*" >&2
    exit 1
}

printHelpAndExit() {
    echo "Usage: ${0##*/} -v VCF -s SAMPLE_INFO"
    echo "-v VCF : path to jointly called VCF (gzipped)"
    echo "-s SAMPLE_INFO : JSON string with sample information"
    exit "$1"
}
while getopts "v:s:" opt; do
    case $opt in
        v) joint_called_vcf="$OPTARG"
           joint_called_vcf_tbi="$OPTARG.tbi"
        ;;
        s) sample_info=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "============================="
echo "Jointly-called VCF: $joint_called_vcf"
echo "Jointly-called VCF index: $joint_called_vcf_tbi"
echo ""
echo "Sample info: $sample_info"
echo "============================="

if [ -z "$joint_called_vcf" ]
then
    echoerr "Annotated VCF missing"
fi

if [ -z "$joint_called_vcf_tbi" ]
then
    echoerr "Annotated VCF index missing"
fi

if [ -z "$sample_info" ]
then
    echoerr "Sample info is missing"
fi



#SCRIPT_LOCATION="/usr/local/bin" # To use in prod
SCRIPT_LOCATION="/Users/alexandervelt/Documents/GitHub/cgap-pipeline-cohort/dockerfiles/regenie/scripts" # To use locally

# Run peddy to infer the ancestry. This will be added to the sample_info json
echo ""
echo "== Run Peddy to infer ancestry =="
sample_info=$(python "$SCRIPT_LOCATION"/run_peddy.py -a "$joint_called_vcf" -s "$sample_info" || exit 1)

# Remove chrM - regenie does not work with it
echo ""
echo "== Removing unsupported chromosomes =="
#bcftools filter "$joint_called_vcf" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY -O z > tmp.no_chrM.vcf.gz || exit 1
bcftools filter "$joint_called_vcf" -r chr1 -O z > tmp.no_chrM.vcf.gz || exit 1
bcftools index -t tmp.no_chrM.vcf.gz || exit 1

# Assign an ID to each variants - existing IDs will be overwritten as these can contain duplicate IDs
echo ""
echo "== Assigning unique ID to each variant =="
#bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' tmp.no_chrM.vcf.gz > tmp.no_chrM.id.vcf || exit 1
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' tmp.no_chrM.vcf.gz -O z > tmp.no_chrM.id.vcf.gz || exit 1
bcftools index -t tmp.no_chrM.id.vcf.gz || exit 1

rm -f tmp.no_chrM.vcf.gz

# Perform variant and sample filtering
echo ""
echo "== Performing variant filtering =="
vcftools --gzvcf tmp.no_chrM.id.vcf.gz \
         --recode \
         --recode-INFO-all \
         --max-missing 0.9 \
         --min-alleles 2 \
         --max-alleles 2 \
         --minQ 90 \
         --minDP 10 \
         --mac 1 --stdout | gzip -c > tmp.no_chrM.id.filtered.recode.vcf.gz || exit 1


echo ""
echo "== Perform Hardy-Weinberg filtering by population =="
python "$SCRIPT_LOCATION"/create_hwe_popmap.py -s "$sample_info" -o tmp.popmap.txt || exit 1
"$SCRIPT_LOCATION"/filter_hwe_by_pop.pl -v tmp.no_chrM.id.filtered.recode.vcf.gz -p tmp.popmap.txt -o tmp.no_chrM.id.hwe.vcf.gz || exit 1
rm -f tmp.no_chrM.id.filtered.recode.vcf.gz

echo ""
echo "== Apply GATK best practice filter =="
python "$SCRIPT_LOCATION"/apply_gatk_filter.py -a tmp.no_chrM.id.hwe.vcf.gz -o joint_called_vcf_filtered.vcf.gz || exit 1

echo ""
echo "== DONE =="

