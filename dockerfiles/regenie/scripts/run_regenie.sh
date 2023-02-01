#!/bin/bash
shopt -s extglob

echoerr() { 
    printf "%s\n" "$*" >&2
    exit 1
}

printHelpAndExit() {
    echo "Usage: ${0##*/} -v VCF -c CASES -e EXCLUDED_GENES -g GENE_ANNOTATIONS -b VC_TESTS -a AAF_BIN"
    echo "-v VCF : path to VEP annotated VCF (gzipped)"
    echo "-c CASES : comma separated list of sample IDs that are affected"
    echo "-a AAF_BIN : specifies the AAF upper bound used to generate burden masks"
    echo "-b VC_TESTS : gene-based tests to use"
    echo "-e EXCLUDED_GENES : comma separated list of genes to exclude from the analysis (no spaces)"
    echo "-g GENE_ANNOTATIONS : gene annotation file from portal"
    exit "$1"
}
while getopts "v:c:g:a:b:e:" opt; do
    case $opt in
        v) annotated_vcf="$OPTARG"
           annotated_vcf_tbi="$OPTARG.tbi"
        ;;
        c) cases=$OPTARG;;
        g) gene_annotations=$OPTARG;;
        a) aaf_bin=$OPTARG;;
        b) vc_tests=$OPTARG;;
        e) excluded_genes=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "============================="
echo "Annotated VCF: $annotated_vcf"
echo "Annotated VCF index: $annotated_vcf_tbi"
echo "Gene annotation file: $gene_annotations"
echo ""
echo "Cases: $cases" #MSA_singleton_99_40-WES,MSA_singleton_782532_01_C_1_S8-WES,MSA_singleton_551913_P1_1_EXO_551913_01-WES,MSA_singleton_541376_01_C_1_S7-WES,MSA_singleton_300683_01_C_1_S6-WES,MSA_singleton_22699_MSA_A_1_S23-WES,MSA_singleton_22692_MSA_A_2_S24-WES,MSA_singleton_22655_MSA_A_2_S21-WES,MSA_singleton_22571_MSA_A_1_S20-WES,MSA_singleton_22561_MSA_A_1_S19-WES,MSA_singleton_22549_MSA_A_2_S18-WES,MSA_singleton_22527_MSA_A_1_S17-WES,MSA_singleton_22514_MSA_A_1_S16-WES,MSA_singleton_22494_MSA_A_1_S15-WES,MSA_singleton_22365_MSA_A_1_S14-WES,MSA_singleton_198677_01_C_1_S5-WES,MSA_singleton_171099_01_MSA_EXO-WES,MSA_singleton_14-49_A_2_S13-WES,MSA_singleton_126939_01_C_1_S2-WES,MSA_singleton_12_18-WES,MSA_singleton_11_46-WES,MSA_singleton_100584_01_C_2_S1-WES,MSA_singleton_07_36-WES,MSA_singleton_07_03-WES,MSA_singleton_04_56-WES,MSA_singleton_04_51-WES,MSA_singleton_03_55-WES,MSA_994562_P1_1_EXO_994562_01-WES,MSA_994562_M1_2_EXO_994562_02-WES,MSA_994562_F1_1_EXO_994562_03-WES,MSA_985648_S1_EXO_S1_MSA-WES,MSA_985648_P1_EXO_Proband_MSA-WES,MSA_985648_M_EXO_M_MSA-WES,MSA_828071_S3_2_EXO_828071_06-WES,MSA_828071_S1_1_EXO_828071_04-WES,MSA_828071_P1_2_EXO_828071_01-WES
echo ""
echo "Genes exluded from analysis: $excluded_genes"
echo "============================="

if [ -z "$annotated_vcf" ]
then
    echoerr "Annotated VCF missing"
fi

if [ -z "$annotated_vcf_tbi" ]
then
    echoerr "Annotated VCF index missing"
fi

if [ -z "$cases" ]
then
    echoerr "Cases missing"
fi

if [ -z "$gene_annotations" ]
then
    echoerr "Gene annotations missing"
fi

if [ -z "$aaf_bin" ]
then
    echoerr "AAF_BIN missing"
fi

if [ -z "$vc_tests" ]
then
    echoerr "Gene-based tests missing. The following are supported: burden,skat,skato,skato-acat,acatv,acato,acato-full."
fi

if [ "$vc_tests" = "burden" ]
then
    vc_tests=""
fi

SCRIPT_LOCATION="/usr/local/bin" # To use in prod
#SCRIPT_LOCATION="/Users/alexandervelt/Documents/GitHub/cgap-pipeline-cohort/dockerfiles/regenie/scripts" # To use in prod

# Remove chrM - regenie does not work with it
echo ""
echo "== Removing unsupported chromosomes =="
bcftools filter "$annotated_vcf" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY -O z > tmp.no_chrM.vcf.gz || exit 1
#bcftools filter "$annotated_vcf" -r chr1 -O z > tmp.no_chrM.vcf.gz || exit 1
bcftools index -t tmp.no_chrM.vcf.gz || exit 1

# Assign an ID to each variants - existing IDs will be overwritten as these can contain duplicate IDs
echo ""
echo "== Assigning unique ID to each variant =="
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' tmp.no_chrM.vcf.gz > tmp.no_chrM.id.vcf || exit 1

# Perform variant and sample filtering
echo ""
echo "== Perform variant filtering =="
vcftools --vcf tmp.no_chrM.id.vcf \
         --recode \
         --recode-INFO-all \
         --out tmp.no_chrM.id.filtered \
         --hwe 0.001 \
         --max-missing 0.9 \
         --min-alleles 2 \
         --max-alleles 2 \
         --minQ 90 \
         --minDP 10 \
         --mac 1 || exit 1

echo ""
echo "== Apply GATK best practice filter =="
python "$SCRIPT_LOCATION"/apply_gatk_filter.py -a tmp.no_chrM.id.filtered.recode.vcf -o tmp.no_chrM.id.filtered.recode.gatk.vcf || exit 1

bgzip -c tmp.no_chrM.id.filtered.recode.gatk.vcf > tmp.no_chrM.id.filtered.recode.gatk.vcf.gz || exit 1
tabix -p vcf tmp.no_chrM.id.filtered.recode.gatk.vcf.gz || exit 1
mv tmp.no_chrM.id.filtered.recode.gatk.vcf regenie_input_source.vcf


# Create BGEN for input to regenie
echo ""
echo "== Create BGEN file with index =="
plink2 --export bgen-1.2 'bits=8' --out regenie_input --vcf tmp.no_chrM.id.filtered.recode.gatk.vcf.gz || exit 1

# Create BGEN index
bgenix -g regenie_input.bgen -index -clobber || exit 1

# Remove temporary files
rm -f tmp*

# Create the phenotype file
# This will create the file 'regenie_input.phenotype'. The sample file is created together with the bgen file.
echo ""
echo "== Create phenotype file =="
python "$SCRIPT_LOCATION"/create_phenotype.py -s regenie_input.sample -o regenie_input.phenotype -c "$cases" || exit 1

# Create files for gene-level testing
# This will create the files 'regenie_input.annotation', 'regenie_input.set_list', 'regenie_input.masks'
echo ""
echo "== Create mask files =="
python "$SCRIPT_LOCATION"/create_mask_files.py -a regenie_input_source.vcf || exit 1

echo ""
echo "== Regenie step 1 =="
regenie --step 1 \
        --bgen regenie_input.bgen \
        --sample regenie_input.sample \
        --phenoFile regenie_input.phenotype \
        --bsize 100 \
        --bt \
        --out regenie_result_step1 || exit 1


echo ""
echo "== Regenie Step 2 - Variant statistics =="
regenie --step 2 \
        --bgen regenie_input.bgen \
        --sample regenie_input.sample \
        --phenoFile regenie_input.phenotype \
        --bsize 200 \
        --bt \
        --firth --approx \
        --pred regenie_result_step1_pred.list \
        --out regenie_result_step2_variant || exit 1

echo ""
echo "== Create variant level result file and Higlass VCF =="

python "$SCRIPT_LOCATION"/create_variant_result_file.py -r regenie_result_step2_variant_Y1.regenie \
                                      -a regenie_input_source.vcf \
                                      -c "$cases" \
                                      -o variant_level_results.txt \
                                      -e higlass_variant_tests.vcf || exit 1

echo ""
echo "== Create multilevel version of the Higlass VCF =="
# This needs "pip install cgap-higlass-data"
cat higlass_variant_tests.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_variant_tests.sorted.vcf
create-cohort-vcf -i higlass_variant_tests.sorted.vcf \
                  -o higlass_variant_tests.multires.vcf \
                  -c fisher_ml10p_control \
                  -q True || exit 1
cat higlass_variant_tests.multires.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_variant_tests.multires.sorted.vcf || exit 1
bgzip -c higlass_variant_tests.multires.sorted.vcf > higlass_variant_tests.multires.vcf.gz || exit 1
tabix -p vcf higlass_variant_tests.multires.vcf.gz || exit 1

echo ""
echo "== Regenie Step 2 - Gene-level statistics =="

regenie --step 2 \
        --bgen regenie_input.bgen \
        --sample regenie_input.sample \
        --phenoFile regenie_input.phenotype \
        --anno-file regenie_input.annotation \
        --set-list regenie_input.set_list \
        --mask-def regenie_input.masks \
        --bsize 200 \
        --bt \
        --firth --approx \
        --minMAC 0.5 \
        --pred regenie_result_step1_pred.list \
        --check-burden-files \
        --strict-check-burden \
        --verbose \
        --aaf-bins "$aaf_bin" \
        --vc-maxAAF "$aaf_bin" \
        --exclude-setlist "$excluded_genes" \
        --out regenie_result_step2_gene \
        --write-mask-snplist \
        --write-mask \
        --vc-tests "$vc_tests" || exit 1



echo ""
echo "== Create gene level Higlass file =="

# Extract gene annotations if necessary
if file --mime-type "$gene_annotations" | grep -q gzip$; then
  gunzip -c "$gene_annotations" > gene_annotations.tsv
else
  cp "$gene_annotations" gene_annotations.tsv
fi

python "$SCRIPT_LOCATION"/create_higlass_gene_file.py -r regenie_result_step2_gene_Y1.regenie \
                                   -g gene_annotations.tsv \
                                   -s regenie_result_step2_gene_masks.snplist \
                                   -a "$aaf_bin" \
                                   -o higlass_gene_tests.vcf || exit 1

cat higlass_gene_tests.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > higlass_gene_tests.sorted.vcf || exit 1
bgzip -c higlass_gene_tests.sorted.vcf > higlass_gene_tests.sorted.vcf.gz || exit 1
tabix -p vcf higlass_gene_tests.sorted.vcf.gz || exit 1


echo ""
echo "== DONE =="

