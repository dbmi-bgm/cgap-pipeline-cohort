cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: joint_called_vcf
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the jointly called vcf.gz file

  - id: reference
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: expect the path to the fa reference file

  - id: sample_info
    type: string
    doc: encoded JSON string containing sample information

outputs:
  joint_called_vcf_filtered:
    type: File
    outputSource: filtering/joint_called_vcf_filtered

steps:
  split_multiallelics:
    run: split_multiallelics.cwl
    in:
      input:
        source: joint_called_vcf
      reference:
        source: reference
    out: [output]

  filtering:
    run: filtering.cwl
    in:
      joint_called_vcf:
        source: split_multiallelics/output
      sample_info:
        source: sample_info
    out: [joint_called_vcf_filtered]

doc: |
  run run_filtering.sh to filter the jointly-called VCF
