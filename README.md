# meta-METAL_nextflow_dsl2
A nextflow DSL2 based pipeline to run METAL based meta-analysis and subsequent plotting of results.

The inspiration, some of the workflow logic and parts of the Dockerfile have been copied and adapted from [Lifebit-ai/metagwas](https://github.com/lifebit-ai/metagwas/blob/stable/Dockerfile) repository.
In contrast to the original Lifebit repository, this workflow uses DSL2 code to run.

To build the docker container locally and run R inside it for testing purposes:

```
docker build -t meta_metal .
docker run -ti --rm meta_metal

## For pushing docker image to Dockerhub
docker tag meta_metal wuennemannflorian/meta_metal:latest
docker push wuennemannflorian/meta_metal:latest
```

Building a singularity sif file from the docker container on Dockerhub:

```
## From local Docker build
singularity build meta_metal.2021_07.sif docker-daemon://meta_metal:latest

## From public Dockerhub
singularity pull -F meta_metal.2021_07.sif docker://wuennemannflorian/meta_metal:latest

## For VEP singularity container
singularity pull -F vep_140.3.sif docker://ensemblorg/ensembl-vep:release_104.3

## Annovar
singularity pull -F annovar.sif docker://bioinfochrustrasbourg/annovar:2018Apr16
```