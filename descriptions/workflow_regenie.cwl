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

  - id: gene_annotations
    type: File
    doc: gene annotations from portal

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

  - id: af_threshold_higlass
    type: float
    default: 0.03
    doc: AF threshold for rare variants. Influences the Higlass result file

  - id: excluded_genes
    type: string
    doc: genes to exclude

outputs:
  variant_level_results:
    type: File
    outputSource: regenie/variant_level_results

  regenie_gene_results:
    type: File
    outputSource: regenie/regenie_gene_results

  higlass_variant_result:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: regenie/higlass_variant_result
  
  higlass_gene_result:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: regenie/higlass_gene_result

  coverage:
    type: File
    outputSource: regenie/coverage

steps:
  regenie:
    run: regenie.cwl
    in:
      annotated_vcf:
        source: annotated_vcf
      sample_info:
        source: sample_info
      gene_annotations:
        source: gene_annotations
      aaf_bin:
        source: aaf_bin
      vc_tests:
        source: vc_tests
      high_cadd_threshold:
        source: high_cadd_threshold
      af_threshold_higlass:
        source: af_threshold_higlass
      excluded_genes:
        source: excluded_genes
    out: [variant_level_results, regenie_gene_results, higlass_variant_result, higlass_gene_result, coverage]

doc: |
  run run_regenie.sh to create statistical analysis results and Higlass files
