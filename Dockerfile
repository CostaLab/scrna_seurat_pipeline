FROM ubuntu:latest
ARG DEBIAN_FRONTEND=noninteractive

## THIS IS FOR Seurat 4 and related softwares


RUN useradd -m docker && echo "docker:docker" | chpasswd
#&& adduser docker sudo

WORKDIR /app

RUN apt-get update && \
    apt-get install -y git

RUN apt-get -y pandoc

RUN apt-get install -y \
    python3.10 \
    python3-pip

RUN apt install -y python3-rpy2
RUN apt install -y python3-jinja2


RUN pip install --break-system-packages --upgrade pip

RUN pip install --break-system-packages grip
RUN pip install --break-system-packages magic-impute
RUN pip install --break-system-packages leidenalg

RUN apt-get install -y r-base
RUN apt-get install -y libarchive-dev
RUN apt-get install -y libharfbuzz-dev libfribidi-dev #devtools
RUN apt-get install -y libcurl4-openssl-dev   #org.Mm.eg.db
RUN apt-get install -y libssl-dev             #org.Mm.eg.db
RUN apt-get install -y libxml2-dev            #org.Mm.eg.db
RUN apt-get install -y libfontconfig1-dev     #clusterprofiler
RUN apt-get install -y libmagick++-dev        #ComplexHeatmap
RUN apt-get install -y libfftw3-dev           #celda
RUN apt-get install -y libudunits2-dev        #schex
RUN apt-get install -y libsqlite3-dev         #schex
RUN apt-get install -y libgdal-dev            #schex

## install cran R packages
RUN  R -e "install.packages('devtools', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('remotes', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "remotes::install_version('SeuratObject', '4.1.4', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN  R -e "remotes::install_version('Seurat', '4.4.0', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN  R -e "install.packages('optparse', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('futile.logger', dependencies=TRUE, repos='http://cloud.r-project.org/')"
#RUN  R -e "install.packages('Seurat', dependencies=TRUE, repos='http://cloud.r-project.org/', version='4.3.0')"
RUN  R -e "install.packages('dplyr', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('future.apply', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('WriteXLS', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('clustree', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Matrix', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('data.table', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Hmisc', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('foreach', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('doParallel', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('glue', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('openxlsx', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('rmarkdown', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('reshape2', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('circlize', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('kableExtra', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('assertthat', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('mclust', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('shinythemes', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('systemfonts', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('igraph', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('proj4', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Cairo', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggalt', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('urltools', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('downloadthis', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('jsonlite', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('crayon', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('SoupX', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('rcompanion', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggridges', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggsci', dependencies=TRUE, repos='http://cloud.r-project.org/')"
#RUN  R -e "install.packages('Rmagic', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('forcats', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('configr', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('archive', dependencies=TRUE, repos='http://cloud.r-project.org/')"


## install Bioconductor R packages
RUN R -e "BiocManager::install('org.Mm.eg.db', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('clusterProfiler', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('org.Hs.eg.db', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('ComplexHeatmap', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('EnhancedVolcano', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('ReactomePA', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('msigdbr', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('limma', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('celda', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('progeny', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('scran', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('GSEABase', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('Nebulosa', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('SingleCellExperiment', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('schex', dependencies=TRUE, ask=FALSE,update=FALSE)"

## install github R packages, might issues too many calls of api.github.com
RUN R -e "devtools::install_github('mahmoudibrahim/genesorteR', dependencies=TRUE, upgrade=FALSE)"
RUN R -e "devtools::install_github('ggjlab/scMCA', dependencies=TRUE, upgrade=FALSE)"
RUN R -e "devtools::install_github('immunogenomics/harmony@7b5f4aa24c5ab34e7015ac2291821102d42cbe63', dependencies=TRUE, upgrade=FALSE)"
RUN R -e "devtools::install_github('ggjlab/scHCL', dependencies=TRUE, upgrade=FALSE)"
RUN R -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder@0710f5a0abae0fdd227fb43a72aed41e94b25152', dependencies=TRUE, upgrade=FALSE)" ## for seurat4
RUN R -e "devtools::install_github('satijalab/seurat-wrappers@796aa878d4ae98776c193cf7c62ca806ba016c0f', dependencies=TRUE, upgrade=FALSE)" #for seurat4
RUN R -e "devtools::install_github('immunogenomics/presto', dependencies=TRUE, upgrade=FALSE)"
#RUN R -e "devtools::install_github('cran/Rmagic',dependencies=TRUE, upgrade=FALSE)"
RUN "git clone https://github.com/KrishnaswamyLab/MAGIC.git && cd Rmagic && R CMD INSTALL ." ##3.0 version
