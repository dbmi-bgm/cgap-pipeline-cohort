#!/usr/bin/env cwl-runner

cwlVersion: v1.0
baseCommand: gather_results.sh
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
  af_threshold_higlass:
    type: float
    inputBinding:
      prefix: -r
      position: 5
  regenie_variant_results:
    type: File
    inputBinding:
      prefix: -b
      position: 6
  regenie_gene_results:
    type: File
    inputBinding:
      prefix: -c
      position: 7
  regenie_gene_results_snplist:
    type: File
    inputBinding:
      prefix: -d
      position: 8
outputs:
  variant_level_results:
    type: File
    outputBinding:
      glob: variant_level_results.txt.gz
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
  coverage:
    type: File
    outputBinding:
      glob: coverage.bw

hints:
  - dockerPull: ACCOUNT/cohort_higlass:VERSION
    class: DockerRequirement
class: CommandLineTool
