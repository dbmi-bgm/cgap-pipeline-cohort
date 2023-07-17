#!/bin/bash
shopt -s extglob

echoerr() { 
    printf "%s\n" "$*" >&2
    exit 1
}

printHelpAndExit() {
    echo "Usage: ${0##*/} -v VCF -s SAMPLE_INFO -e EXCLUDED_GENES -b VC_TESTS -a AAF_BIN"
    echo "-v VCF : path to VEP annotated VCF (gzipped)"
    echo "-s SAMPLE_INFO : JSON string with sample information"
    echo "-a AAF_BIN : specifies the AAF upper bound used to generate burden masks"
    echo "-b VC_TESTS : gene-based tests to use"
    echo "-e EXCLUDED_GENES : comma separated list of genes to exclude from the analysis (no spaces)"
    exit "$1"
}
while getopts "v:s:g:a:b:c:r:e:" opt; do
    case $opt in
        v) annotated_vcf="$OPTARG"
           annotated_vcf_tbi="$OPTARG.tbi"
        ;;
        s) sample_info=$OPTARG;;
        a) aaf_bin=$OPTARG;;
        b) vc_tests=$OPTARG;;
        c) high_cadd_threshold=$OPTARG;;
        e) excluded_genes=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "============================="
echo "Annotated VCF: $annotated_vcf"
echo "Annotated VCF index: $annotated_vcf_tbi"
echo ""
echo "Sample info: $sample_info"
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

if [ -z "$sample_info" ]
then
    echoerr "Sample info is missing"
fi

if [ -z "$aaf_bin" ]
then
    echoerr "AAF_BIN missing"
fi

if [ -z "$vc_tests" ]
then
    echoerr "Gene-based tests missing. The following are supported: burden,skat,skato,skato-acat,acatv,acato,acato-full."
fi

if [ -z "$high_cadd_threshold" ]
then
    echoerr "High CADD threshold missing."
fi

if [ "$vc_tests" = "burden" ]
then
    vc_tests=""
fi

SCRIPT_LOCATION="/usr/local/bin" # To use in prod
#SCRIPT_LOCATION="/Users/alexandervelt/Documents/GitHub/cgap-pipeline-cohort/dockerfiles/regenie/scripts" # To use locally


# Create BGEN for input to regenie
echo ""
echo "== Create BGEN file with index =="
plink2 --export bgen-1.2 'bits=8' --out regenie_input --vcf "$annotated_vcf" || exit 1

# Create BGEN index
bgenix -g regenie_input.bgen -index -clobber || exit 1

# Remove temporary files
rm -f tmp*

# Create the phenotype file
# This will create the file 'regenie_input.phenotype'. The sample file is created together with the bgen file.
echo ""
echo "== Create phenotype file =="
python "$SCRIPT_LOCATION"/create_phenotype.py -s regenie_input.sample -o regenie_input.phenotype -c "$sample_info" || exit 1

# Create files for gene-level testing
# This will create the files 'regenie_input.annotation', 'regenie_input.set_list', 'regenie_input.masks'
echo ""
echo "== Create mask files =="
python "$SCRIPT_LOCATION"/create_mask_files.py -a "$annotated_vcf" -c "$high_cadd_threshold" || exit 1


echo ""
echo "== Regenie step 1 =="

# Exract at most 500k high-quality variants for step 1
plink2 --bgen regenie_input.bgen 'ref-last' --sample regenie_input.sample --maf 0.01 --mac 10 --geno 0.1 --mind 0.1 --out qc_pass --snps-only --export bgen-1.2 'bits=8' || exit 1
plink2 --bgen qc_pass.bgen 'ref-last' --sample qc_pass.sample --thin-count 500000 --write-snplist --write-samples --no-id-header --out qc_pass_500k || exit 1


regenie --step 1 \
        --bgen regenie_input.bgen \
        --extract qc_pass_500k.snplist \
        --keep qc_pass_500k.id \
        --sample regenie_input.sample \
        --phenoFile regenie_input.phenotype \
        --bsize 100 \
        --bt \
        --lowmem \
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
echo "== Compress result files =="

gzip -c regenie_result_step2_variant_Y1.regenie > regenie_result_variant_Y1.txt.gz
gzip -c regenie_result_step2_gene_Y1.regenie > regenie_result_gene_Y1.txt.gz
gzip -c regenie_result_step2_gene_masks.snplist > regenie_result_gene_masks.snplist.gz


echo ""
echo "== DONE =="

