#!/bin/bash

### Job name
#SBATCH -J viz_ex
#SBATCH -o logs/%j_job.log
#SBATCH -e logs/%j_job.log

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 24:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=48G

### OpenMP threads
#SBATCH --cpus-per-task=8

################################################################
module load R
module load scRNA
################################################################

# could also define it here instead of taking it as an arg
proj_name="$1"
# data dir (where your results were saved)
# the way it is set up below -o for output will also save your report
# on that same directory path
data_path="/data/EXAMPLE/exp/scRNA/some_project/scrna_seurat_pipeline_results"

date
bash ./viz/create_report.sh \
	-a "T.Maié" \
	-p ${proj_name} \
	-m \
	-s "${data_path}/save_${proj_name}" \
	-c "./conf/config_${proj_name}.R" \
	-o "${data_path}/report_${proj_name}" \
	-ch "${data_path}/charts_${proj_name}" \
	-e ./external/Human_and_mouse_cell_markers-Markers.tsv
date
