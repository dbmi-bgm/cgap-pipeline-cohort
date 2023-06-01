cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_vcf
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the vcf gz file

  - id: reference
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: expect the path to the fa reference file

  - id: regions
    type: File
    doc: expect the path to the file defining regions

  - id: vep
    type: File
    secondaryFiles:
      - ^^^.plugins.tar.gz
    doc: expect the path to the vep tar gz

  - id: dbnsfp
    type: File
    secondaryFiles:
      - .tbi
      - ^.readme.txt
    doc: expect the path to the vcf gz file with dbNSFP annotations

  - id: spliceai_snv
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the vcf gz file with SpliceAI for SNVs

  - id: spliceai_indel
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the vcf gz file with SpliceAI for indels

  - id: gnomad
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the vcf gz file with gnomAD annotations

  - id: gnomad2
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the vcf gz file with gnomAD v2.1 exome annotations

  - id: CADD_snv
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the tsv gz file with CADD Phred annotations for SNVs

  - id: CADD_indel
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the tsv gz file with CADD Phred annotations for indels

  - id: nthreads
    type: int
    doc: number of threads used to run parallel

  - id: version
    type: string
    default: "101"
    doc: vep datasource version

  - id: assembly
    type: string
    default: "GRCh38"
    doc: genome assembly version

outputs:
  annotated_vcf:
    type: File
    outputSource: vep_annot/output

steps:
  bcftools_norm_multiallelics:
    run: bcftools_norm_multiallelics.cwl
    in:
      input:
        source: input_vcf
      reference:
        source: reference
    out: [output]

  vep_annot:
    run: vep_annot.cwl
    in:
      input:
        source: bcftools_norm_multiallelics/output
      reference:
        source: reference
      regions:
        source: regions
      vep:
        source: vep
      dbnsfp:
        source: dbnsfp
      spliceai_snv:
        source: spliceai_snv
      spliceai_indel:
        source: spliceai_indel
      gnomad:
        source: gnomad
      gnomad2:
        source: gnomad2
      CADD_snv:
        source: CADD_snv
      CADD_indel:
        source: CADD_indel
      nthreads:
        source: nthreads
      version:
        source: version
      assembly:
        source: assembly
    out: [output]


doc: |
  run bcftools_norm_multiallelics |
  run vep
