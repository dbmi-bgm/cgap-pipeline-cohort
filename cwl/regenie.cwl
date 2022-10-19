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
  cases:
    type: string
    inputBinding:
      prefix: -c
      position: 2
  gene_annotations:
    type: File
    inputBinding:
      prefix: -g
      position: 3
  aaf_bin:
    type: int
    inputBinding:
      prefix: -a
      position: 4
  vc_tests:
    type: string
    inputBinding:
      prefix: -b
      position: 5
  excluded_genes:
    type: string
    inputBinding:
      prefix: -e
      position: 6
outputs:
  regenie_variant_result:
    type: File
    outputBinding:
      glob: regenie_result_step2_variant_Y1.regenie
  regenie_gene_result:
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

hints:
  - dockerPull: aveit/cgap-regenie:0.1.0
    class: DockerRequirement
class: CommandLineTool
