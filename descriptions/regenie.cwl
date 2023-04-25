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
  gene_annotations:
    type: File
    inputBinding:
      prefix: -g
      position: 3
  aaf_bin:
    type: float
    inputBinding:
      prefix: -a
      position: 4
  vc_tests:
    type: string
    inputBinding:
      prefix: -b
      position: 5
  high_cadd_threshold:
    type: float
    inputBinding:
      prefix: -c
      position: 6
  af_threshold_higlass:
    type: float
    inputBinding:
      prefix: -r
      position: 7
  excluded_genes:
    type: string
    inputBinding:
      prefix: -e
      position: 8
outputs:
  variant_level_results:
    type: File
    outputBinding:
      glob: variant_level_results.txt
  regenie_gene_results:
    type: File
    outputBinding:
      glob: regenie_result_step2_gene_Y1.regenie
  higlass_variant_result:
    type: File
    outputBinding:
      glob: higlass_variant_tests.multires.vcf.gz
    secondaryFiles:
      - .tbi
  higlass_gene_result:
    type: File
    outputBinding:
      glob: higlass_gene_tests.sorted.vcf.gz
    secondaryFiles:
      - .tbi
  annotated_vcf_filtered:
    type: File
    outputBinding:
      glob:   annotated_vcf_filtered.vcf.gz
    secondaryFiles:
      - .tbi
  coverage:
    type: File
    outputBinding:
      glob: coverage.bw

hints:
  - dockerPull: ACCOUNT/regenie:VERSION #aveit/cgap-regenie:0.1.1
    class: DockerRequirement
class: CommandLineTool
