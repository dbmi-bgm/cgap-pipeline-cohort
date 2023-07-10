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

  - id: af_threshold_higlass
    type: float
    default: 0.03
    doc: AF threshold for rare variants. Influences the Higlass result file

  - id: regenie_variant_results
    type: File
    doc: Regenie variant result file

  - id: regenie_gene_results
    type: File
    doc: Regenie gene result file

  - id: regenie_gene_results_snplist
    type: File
    doc: Regenie gene result snplist file

outputs:
  variant_level_results:
    type: File
    outputSource: higlass/variant_level_results

  higlass_variant_result:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: higlass/higlass_variant_result
  
  higlass_gene_result:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: higlass/higlass_gene_result

  coverage:
    type: File
    outputSource: higlass/coverage

steps:
  higlass:
    run: higlass.cwl
    in:
      annotated_vcf:
        source: annotated_vcf
      sample_info:
        source: sample_info
      gene_annotations:
        source: gene_annotations
      aaf_bin:
        source: aaf_bin
      af_threshold_higlass:
        source: af_threshold_higlass
      regenie_variant_results:
        source: regenie_variant_results
      regenie_gene_results:
        source: regenie_gene_results
      regenie_gene_results_snplist:
        source: regenie_gene_results_snplist

    out: [variant_level_results, higlass_variant_result, higlass_gene_result, coverage]

doc: |
  Create all the result files from the analysis
