#!/bin/bash

### Job name
#SBATCH -J scrna_test_2026-02-25_10:23:08.834907

### Set logs, remember to create logs dir
#SBATCH -e logs/scrna_test_2026-02-25_10:23:08.834907.txt
#SBATCH -o logs/scrna_test_2026-02-25_10:23:08.834907.txt

### Time to execute, e. g. 15 min 30 sec
#SBATCH -t 96:00:00

### Job memory needs per node, e. g. 1 GB
#SBATCH --mem=180G

### OpenMP threads
#SBATCH --cpus-per-task=24

################################################################
# PATH
#if [ -r /usr/local_host/etc/bashrc ]; then
#   . /usr/local_host/etc/bashrc
#fi

#export PATH=/usr/bin:$PATH
#export PATH=/usr/local_host/bin:$PATH
################################################################
#module load scRNA
#source ~/miniconda3/bin/activate
#conda activate Seurat5
################################################################

# could also define it here instead of taking it as an arg
proj_name=$1
# data dir (where your results will be saved)

data_path="/data/EXAMPLE/exp/scRNA/some_project/scrna_seurat_pipeline_results"
data_path=`pwd`

date
## 50 cores run, future memory
mkdir -p ${data_path}/${proj_name}
ln -s ${data_path}/conf/config_${proj_name}.R ${data_path}/${proj_name}
Rscript data_factory.R \
  -n 24 \
  --MaxMemMega=180000 \
  -z "lz4" \
  -c "./conf/config_${proj_name}.R" \
  -s "${data_path}/${proj_name}/save" \
  -e "${data_path}/${proj_name}/charts" \
  -a seurat_clusters \
  --allinone=TRUE \
  --offline_all_databases=TRUE \
  --nFeatureRNAfloor=400 \
  --nCountRNAfloor=0 \
  --nCountRNAceiling=40000 \
  --pct_mitoceiling=20 \
  --pct_riboceiling=100 \
  --countmatrixformat='10X' \
  --dims4Integrate="1:30" \
  --Dims4FindNeighbors="1:50" \
  --harmony_dim="1:50" \
  -r 0.5

date
