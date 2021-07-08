FROM r-base:4.1.0

## Install 
RUN R -e "install.packages(c('qqman','data.table'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

##  Download and unpack METAL
RUN wget http://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz \
    && tar -xzvf Linux-metal.tar.gz \
    && rm Linux-metal.tar.gz

## Add METAL executable to PATH
ENV PATH /generic-metal:$PATH

## Set working directory to /data
WORKDIR /data/