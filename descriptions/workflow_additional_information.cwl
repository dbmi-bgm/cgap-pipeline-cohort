cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: annotated_vcf
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the jointly called, filtered vcf.gz file

  - id: sample_info
    type: string
    doc: encoded JSON string containing sample information

outputs:
  
  variant_details:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: additional_information/variant_details
  
steps:
  additional_information:
    run: additional_information.cwl
    in:
      annotated_vcf:
        source: annotated_vcf
      sample_info:
        source: sample_info
    out: [variant_details]

doc: |
  Creates a file with additional information for each variant of the cohort
