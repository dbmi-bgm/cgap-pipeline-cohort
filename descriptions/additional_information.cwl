#!/usr/bin/env cwl-runner

cwlVersion: v1.0
baseCommand: create_variant_details_file.sh
requirements:
  InlineJavascriptRequirement: {}
inputs:
  annotated_vcf:
    type: File
    inputBinding:
      prefix: -a
      position: 1
    secondaryFiles:
      - .tbi
  sample_info:
    type: string
    inputBinding:
      prefix: -s
      position: 2
outputs:
  variant_details:
    type: File
    outputBinding:
      glob: variant_details.vcf.gz
    secondaryFiles:
      - .tbi

hints:
  - dockerPull: ACCOUNT/cohort_higlass:VERSION
    class: DockerRequirement
class: CommandLineTool
