FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ARG GITHUB_PAT
## THIS IS FOR Seurat 5 and related softwares


RUN useradd -m docker && echo "docker:docker" | chpasswd
#&& adduser docker sudo

WORKDIR /app


RUN apt-get update && \
    apt-get install -y git

RUN apt-get install -y curl

RUN apt-get install -y pandoc

RUN apt-get install -y \
    python3.10 \
    python3.10-dev 



RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip


RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 2
#update-alternatives --config python3


RUN python3 -m pip install jinja2==3.1.6
RUN pip install  grip==4.6.2
RUN pip install  magic-impute==3.0.0
RUN pip install  leidenalg==0.11.0
RUN pip install --force-reinstall  pandas  ## To prevent load pandas failure
RUN python3.11 -c "import pandas as pd; print('pandas OK:', pd.__version__)"



RUN  curl https://cloud.r-project.org/src/base/R-4/R-4.4.0.tar.gz > /tmp/R-4.4.0.tar.gz
RUN cd /tmp && tar -xzf R-4.4.0.tar.gz
WORKDIR /tmp/R-4.4.0


RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
     \
    tzdata \
    libreadline-dev \
    libx11-dev \
    libxt-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libpcre2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    libcairo2-dev \
    && rm -rf /var/lib/apt/lists/*


RUN ./configure \
    --prefix=/opt/R/4.4.0 \
    --enable-R-shlib \
    --with-blas \
    --with-lapack \
    && make -j$(nproc) \
    && make install

ENV PATH="/opt/R/4.4.0/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/R/4.4.0/lib/R/lib:${LD_LIBRARY_PATH:-}"

RUN R --version



RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpcre2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade pip setuptools wheel

RUN python3 -m pip install rpy2==3.6.4

RUN apt-get update && apt-get install -y libarchive13 libarchive-tools
RUN apt-get install -y libharfbuzz-dev #devtools
RUN apt-get install -y libfribidi-dev #devtools
RUN apt-get install -y libcurl4-openssl-dev   #org.Mm.eg.db
RUN apt-get install -y libssl-dev             #org.Mm.eg.db
RUN apt-get install -y libxml2-dev            #org.Mm.eg.db
RUN apt-get install -y libfontconfig1-dev     #clusterprofiler
RUN apt-get install -y libmagick++-dev        #ComplexHeatmap
RUN apt-get install -y libfftw3-dev           #celda
RUN apt-get install -y libudunits2-dev        #schex
RUN apt-get install -y libsqlite-dev         #schex
RUN apt-get install -y libgdal-dev            #schex
RUN apt-get install -y cmake
RUN apt-get install -y libgmp-dev

RUN apt-get install -y libarchive-dev
RUN  R -e "install.packages('archive', version='1.1.12', dependencies=TRUE, repos='http://cloud.r-project.org/')"



## install cran R packages
RUN  R -e "install.packages('devtools', version='2.4.6', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('remotes', version='2.5.0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('usethis', version='3.2.1', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "remotes::install_version('SeuratObject', version='5.3.0', \
                                        repos=c(\
                                       'https://satijalab.r-universe.dev',\
                                       'https://cloud.r-project.org'\
                                          )\
                                    )"
RUN R -e "remotes::install_version('Seurat', version='5.4.0', \
                                        repos=c( \
                                       'https://satijalab.r-universe.dev',\
                                       'https://cloud.r-project.org'\
                                         ) \
                                    )"
RUN  R -e "install.packages('optparse', version='1.7.5', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('futile.logger', version='1.4.9', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('dplyr', version='1.2.0', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('future.apply', version='1.20.1', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('WriteXLS', version='6.8.0', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('clustree', version='0.5.1',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Matrix', version='1.7-4', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('data.table', version='1.17.8', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggplot2', version='3.5.2', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Hmisc', version='5.2-5',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('foreach', version='1.5.2',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('doParallel', version='1.0.17', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('glue', version='1.8.0', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('openxlsx', version='4.2.8.1',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('rmarkdown', version='2.30', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('reshape2', version='1.4.5', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('circlize', version='0.4.17',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('BiocManager', version='1.29.1', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('kableExtra', version='1.4.0', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('assertthat', version='0.2.1', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('mclust', version='6.1.2', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('shinythemes', version='1.2.0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('systemfonts', version='1.3.1',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('igraph', version='2.1.4', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('proj4', version='1.0-14',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('Cairo', version='1.7-0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggalt', version='0.4.0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('urltools', version='1.7.3.1', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('downloadthis', version='0.5.0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('jsonlite', version='2.0.0', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('crayon', version='1.5.3', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('SoupX', version='1.6.2', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('rcompanion', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggridges', version='0.5.7', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('ggsci', version='4.2.0',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('forcats', version='1.0.1',dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('configr', version='0.3.5', dependencies=TRUE, repos='http://cloud.r-project.org/')"
RUN  R -e "install.packages('archive', version='1.1.12', dependencies=TRUE, repos='http://cloud.r-project.org/')"



## install Bioconductor R packages
RUN R -e "BiocManager::install(version = '3.19', ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install('org.Mm.eg.db', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"

RUN R -e "BiocManager::install('org.Hs.eg.db', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('ComplexHeatmap', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('EnhancedVolcano', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('ReactomePA', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('msigdbr', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('limma', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('celda', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('progeny', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('scran', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('GSEABase', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('Nebulosa',version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('SingleCellExperiment', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"
RUN R -e "BiocManager::install('schex', version='3.19',dependencies=TRUE, ask=FALSE,update=FALSE)"

## install github R packages, might issues too many calls of api.github.com
RUN export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('mahmoudibrahim/genesorteR', dependencies=TRUE, upgrade=FALSE)"
RUN export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('ggjlab/scMCA', version='0.2.0',dependencies=TRUE, upgrade=FALSE)"
RUN export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('immunogenomics/harmony', version='1.2.4',dependencies=TRUE, upgrade=FALSE)"
RUN export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('ggjlab/scHCL', version='0.1.1',dependencies=TRUE, upgrade=FALSE)"
RUN export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', version='2.0.6', dependencies=TRUE, upgrade=FALSE)" ## for seurat4
RUN  export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('satijalab/seurat-wrappers', version='0.4.0',dependencies=TRUE, upgrade=FALSE)" #for seurat4
RUN  export GITHUB_PAT=${GITHUB_PAT} &&  R -e "devtools::install_github('immunogenomics/presto', dependencies=TRUE, upgrade=FALSE)"
RUN  export GITHUB_PAT=${GITHUB_PAT} && R -e "devtools::install_github('cran/Rmagic', version='2.0.3',dependencies=TRUE, upgrade=FALSE)"
### -----
RUN R -e "remotes::install_github('ctlab/fgsea')"
RUN R -e "remotes::install_github('YuLab-SMU/tidytree')"
RUN R -e "remotes::install_github('YuLab-SMU/treeio')"
RUN R -e "remotes::install_github('YuLab-SMU/ggtree')"
RUN R -e "BiocManager::install('clusterProfiler', version='3.19', dependencies=TRUE, ask=FALSE,update=FALSE)"

RUN apt-get install -y wget curl
RUN wget "https://bioconductor.org/packages/3.19/bioc/src/contrib/ReactomePA_1.48.0.tar.gz"
RUN R CMD INSTALL ReactomePA_1.48.0.tar.gz
## cache_dir <- tools::R_user_dir(package = "msigdbr", which = "cache") -> /root/.cache/R/msigdbr
#https://zenodo.org/records/15800824/files/msigdb.2025.1.Mm.rds?download=1
#https://zenodo.org/records/15800824/files/msigdb.2025.1.Hs.rds?download=1
RUN cache_dir=/root/.cache/R/msigdbr

RUN curl -4 -L \
     --retry 10 \
     --retry-delay 5 \
     --retry-connrefused \
     -o $cache_dir/msigdb.2025.1.Hs.rds \
     "https://zenodo.org/records/15800824/files/msigdb.2025.1.Hs.rds"

RUN curl -4 -L \
     --retry 10 \
     --retry-delay 5 \
     --retry-connrefused \
     -o $cache_dir/msigdb.2025.1.Mm.rds \
     "https://zenodo.org/records/15800824/files/msigdb.2025.1.Mm.rds"


RUN ln -s /usr/bin/python3 /usr/bin/python
