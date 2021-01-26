#!/bin/bash

### Job name
#SBATCH -J conda 
#SBATCH -e ./conda.txt
#SBATCH -o ./conda.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 100:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=80G

#SBATCH -A  rwth0455


################################################################
# PATH
#if [ -r /usr/local_host/etc/bashrc ]; then
#   . /usr/local_host/etc/bashrc
#fi
#
#export PATH=/usr/bin:$PATH
#export PATH=/usr/local_host/bin:$PATH
#source /home/sz753404/for-se/anaconda3/bin/activate
#conda activate r_env4s3
################################################################
# LIBRARYPATH


Rscript data_factory.R -c conf/configCSCQC.R \
                       -s saveCSCQC/ \
                       -e chartsCSCQC/ \
                       --nFeatureRNAfloor 200 \
		       --nCountRNAceiling 12000 \
		       --nCountRNAfloor 45 \
		       --pct_mitoceiling 20 \
		       --pct_riboceiling 40


## 50 cores run, future memory
# Rscript data_factory.R -n 50  --MaxMemMega=100000


## Preferred setup with a project name (i.e. YOURPROJ)
# Rscript data_factory.R -n 8  --MaxMemMega=180000 -c conf/configYOURPROJ.R -s saveYOURPROJ -e chartsYOURPROJ -a seurat_clusters


# set config file to use, default(conf/config.R) 
#Rscript data_factory.R -c conf/another_config.R 


# set count matrix format, 10X and 10X_h5 
#Rscript data_factory.R -f 10X_h5 


# set cluster resolution(0.5) for metadata slot "seurat_clusters"
#Rscript data_factory.R -r 0.5


# set Find anchors(usually cca) for integration(1:20), dims to use
#Rscript data_factory.R -d 1:20 


# set Find neighbours for(usaully pca after integration) clustering(1:12), dims to use
#Rscript data_factory.R -x 1:12


# set using clusters: "merged_clusters". Do DE or GO analysis
#Rscript data_factory.R -a merged_clusters -c conf/config_merged_clusters.R


# set using clusters: "removed_clusters". Do DE or GO analysis
#Rscript data_factory.R -a removed_clusters -c conf/config_remove_clusters.R


# set using clusters: "remove_recluster". Do DE or GO analysis
#Rscript data_factory.R -a remove_recluster -c conf/config_remove_recluster.R



