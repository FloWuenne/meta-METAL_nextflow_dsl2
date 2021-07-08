# meta-METAL_nextflow_dsl2
A nextflow DSL2 based pipeline to run METAL based meta-analysis and subsequent plotting of results.

The inspiration, some of the workflow logic and parts of the Dockerfile have been copied and adapted from [Lifebit-ai/metagwas](https://github.com/lifebit-ai/metagwas/blob/stable/Dockerfile) repository.
In contrast to the original Lifebit repository, this workflow uses DSL2 code to run.

To build the docker container locally and run R inside it for testing purposes:

```
docker build -t meta_metal .
docker run -ti --rm meta_metal
```

Building a singularity sif file from the docker container on Dockerhub:

```
singularity pull meta_metal.2021_07.sif docker://wuennemannflorian/meta_metal:latest
```