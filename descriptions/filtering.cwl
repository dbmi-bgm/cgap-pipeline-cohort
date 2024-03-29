#!/usr/bin/env cwl-runner

cwlVersion: v1.0
baseCommand: run_filtering.sh
requirements:
  InlineJavascriptRequirement: {}
inputs:
  joint_called_vcf:
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
outputs:
  joint_called_vcf_filtered:
    type: File
    outputBinding:
      glob: joint_called_vcf_filtered.vcf.gz
    secondaryFiles:
      - .tbi
hints:
  - dockerPull: ACCOUNT/cohort_filtering:VERSION
    class: DockerRequirement
class: CommandLineTool
