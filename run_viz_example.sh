#!/bin/bash

### Job name
#SBATCH -J viz_ex
#SBATCH -o logs/%j_job.log
#SBATCH -e logs/%j_job.log

### Time your job needs to execute, e. g. 15 min 30 sec #SBATCH -t 24:00:00
### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=48G

### OpenMP threads
#SBATCH --cpus-per-task=8

################################################################
module load R
module load scRNA
#source ~/miniconda3/bin/activate
#conda activate Seurat3
################################################################

# could also define it here instead of taking it as an arg
proj_name="$1"
# data dir (where your results were saved)
# the way it is set up below -o for output will also save your report
# on that same directory path
data_path="/data/EXAMPLE/exp/scRNA/some_project/scrna_seurat_pipeline_results"
#data_path="/home/sz753404/data/test/pipelines/scrna_seurat_pipeline_unify"

#cluster="removed_clusters"
#cluster="remove_recluster"
#cluster="merged_clusters"
#cluster="annotation"
#cluster="singleton"
cluster="seurat_clusters"

FUNCS=(
  #Singleton
 QC
 AmbientRNA
 DEs
 Clusters
 Clusters_harmony
 Clusters_seurat
 EXT_MARKERS
 DEGO
 progeny
 hallmark
 KEGG
 Reactome
 progeny_stage
 DEGO_stage
 hallmark_stage
 reactome_stage
 kegg_stage
 DEGO_1v1
 hallmark_1v1
 reactome_1v1
 kegg_1v1
 intUMAPs
  #GET_DATA
)

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }
join_str=`join_by '", "' ${FUNCS[@]}`
json_exe_list="[\"${join_str}\"]"

date
Rscript ./viz/create_report.R \
	-a "Mingbo" \
	-p ${proj_name} \
	-m TRUE\
	-s "${data_path}/${proj_name}/save" \
	-c "./conf/config_${proj_name}.R" \
	-o "${data_path}/${proj_name}/report" \
	-r "${data_path}/${proj_name}/charts" \
	-e ./external/Human_and_mouse_cell_markers-Markers.tsv \
    -d "$cluster" \
    -j "${json_exe_list}" \
    -i "FALSE"
date
