# UNDER ACTIVE DEVELOPMENT

<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Pipeline for statistical analysis of case/control cohort data

This repository contains components for the CGAP pipeline for statistical analysis of case/control cohort data:

  * CWL workflow descriptions
  * CGAP Portal *Workflow* and *MetaWorkflow* objects
  * CGAP Portal *Software*, *FileFormat*, and *FileReference* objects
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using [*cgap pipeline utils*](https://github.com/dbmi-bgm/cgap-pipeline-utils/) `pipeline_deploy`)

The pipeline starts from jointly calles, VEP annotated `vcf` files and produces various result file in text and `vcf` format. See documentation for details.




# cgap-pipeline-cohort

To build the docker image run
```
docker build -t cgap/cgap-regenie:0.1.0
```
