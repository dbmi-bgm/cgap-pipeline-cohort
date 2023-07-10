#!/usr/bin/env cwl-runner

cwlVersion: v1.0
baseCommand: run_regenie.sh
requirements:
  InlineJavascriptRequirement: {}
inputs:
  annotated_vcf:
    type: File
    inputBinding:
      prefix: -v
      position: 1
    secondaryFiles:
      - .tbi
  sample_info:
    type: string
    inputBinding:
      prefix: -s
      position: 2
  aaf_bin:
    type: float
    inputBinding:
      prefix: -a
      position: 3
  vc_tests:
    type: string
    inputBinding:
      prefix: -b
      position: 4
  high_cadd_threshold:
    type: float
    inputBinding:
      prefix: -c
      position: 5
  excluded_genes:
    type: string
    inputBinding:
      prefix: -e
      position: 6
outputs:
  regenie_variant_results:
    type: File
    outputBinding:
      glob: regenie_result_variant_Y1.txt.gz
  regenie_gene_results:
    type: File
    outputBinding:
      glob: regenie_result_gene_Y1.txt.gz
  regenie_gene_results_snplist:
    type: File
    outputBinding:
      glob: regenie_result_gene_masks.snplist.gz

hints:
  - dockerPull: ACCOUNT/cohort_regenie:VERSION #aveit/cgap-regenie:0.1.1
    class: DockerRequirement
class: CommandLineTool
