cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: annotated_vcf
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the jointly called vcf.gz file

  - id: sample_info
    type: string
    doc: encoded JSON string containing sample information

  - id: aaf_bin
    type: float
    default: 0.01
    doc: Max alternate allele frequence to consider for gene tests

  - id: vc_tests
    type: string
    default: "skato,acato-full"
    doc: comma separated list of VC tests to run

  - id: high_cadd_threshold
    type: float
    default: 20.0
    doc: CADD threshold for deleterious variants

  - id: excluded_genes
    type: string
    doc: genes to exclude

outputs:
  regenie_variant_results:
    type: File
    outputSource: regenie/regenie_variant_results

  regenie_gene_results:
    type: File
    outputSource: regenie/regenie_gene_results

  regenie_gene_results_snplist:
    type: File
    outputSource: regenie/regenie_gene_results_snplist
  
steps:
  regenie:
    run: regenie.cwl
    in:
      annotated_vcf:
        source: annotated_vcf
      sample_info:
        source: sample_info
      aaf_bin:
        source: aaf_bin
      vc_tests:
        source: vc_tests
      high_cadd_threshold:
        source: high_cadd_threshold
      excluded_genes:
        source: excluded_genes
    out: [regenie_variant_results, regenie_gene_results, regenie_gene_results_snplist]

doc: |
  run run_regenie.sh to create statistical analysis results from Regenie
