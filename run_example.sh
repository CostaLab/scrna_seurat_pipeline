#!/bin/bash

### Job name
#SBATCH -J proj_ex

### Set logs, remember to create logs dir
#SBATCH -e logs/log_%j.log
#SBATCH -o logs/log_%j.log

### Time to execute, e. g. 15 min 30 sec
#SBATCH -t 96:00:00

### Job memory needs per node, e. g. 1 GB
#SBATCH --mem=180G

### OpenMP threads
#SBATCH --cpus-per-task=16

################################################################
# PATH
if [ -r /usr/local_host/etc/bashrc ]; then
   . /usr/local_host/etc/bashrc
fi

export PATH=/usr/bin:$PATH
export PATH=/usr/local_host/bin:$PATH
################################################################
module load scRNA
################################################################

# could also define it here instead of taking it as an arg
proj_name="$1"
# data dir (where your results will be saved)

data_path="/data/EXAMPLE/exp/scRNA/some_project/scrna_seurat_pipeline_results"


date
## 50 cores run, future memory
mkdir -p ${data_path}/${proj_name}
Rscript data_factory.R -n 40 \
  --MaxMemMega=180000 \
  -c "./conf/config_${proj_name}.R" \
  -s "${data_path}/${proj_name}/save" \
  -e "${data_path}/${proj_name}/charts" \
  -a seurat_clusters \
  -r 0.5

date
