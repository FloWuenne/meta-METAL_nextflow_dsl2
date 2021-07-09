FROM rocker/r-ubuntu:20.04

#RUN apt-get update && apt-get install -y procps

##  Download and unpack METAL
RUN wget http://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz \
    && tar -xzvf Linux-metal.tar.gz \
    && rm Linux-metal.tar.gz

## Install 
RUN install2.r qqman data.table magrittr tidyr dplyr R.utils CMplot

RUN mkdir /opt/bin
COPY bin/* /opt/bin/

## Add METAL executable to PATH
ENV PATH /generic-metal:/opt/bin/:$PATH

RUN find /opt/bin/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.R" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.sh" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.css" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.Rmd" -exec chmod +x {} \;

## Set working directory to /data
WORKDIR /data/