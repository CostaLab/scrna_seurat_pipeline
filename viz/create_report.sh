#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage:"
    echo "'create_report.sh [PROJECT_NAME] [true/false](make report elements from scratch)'"
    exit 0
fi

RED='\033[0;31m'
NC='\033[0m' # No Color

PROJ=$1
MAKE_ELEMENTS=$2
SAVED_DATA="../save${PROJ}"
EXT_ANNOT="../external/Human_and_mouse_cell_markers-Markers.tsv"
MAKE_SINGLE_FILE="TRUE"

FUNCS=(
  #Singleton
  QC
  DEs
  Clusters
  Clusters_harmony
  Clusters_seurat
  EXT_MARKERS
  DEGO
  hallmark
  KEGG
  Reactome
  # DEGO_stage
  # hallmark_stage
  # reactome_stage
  # kegg_stage
  # DEGO_1v1
  # hallmark_1v1
  # reactome_1v1
  # kegg_1v1
  intUMAPs
  GET_DATA
)



#!!!!!!!!------clusters to choose-------------
# In general, we choose seurat_clusters,
# If you are using removed or merged clusters,
# choose the following:
            # seurat_clusters
            # merged_clusters
            # removed_clusters
            # remove_recluster

#cluster="removed_clusters"
#cluster="remove_recluster"
#cluster="merged_clusters"
#cluster="annotation"

cluster="seurat_clusters"
#cluster="singleton"

#!!!!!!!!!----------------------------------

echo $PROJ
echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p report${PROJ}/data
cp -r ../charts${PROJ}/* report${PROJ}/data
python code_generator.py -c $cluster -cf "conf/config${PROJ}.R" -b "../" -p "${PROJ}"
grip --export report${PROJ}/index.md

if [[ "$MAKE_ELEMENTS" == true ]]; then
  Rscript make_report_elements.R --proj=$PROJ --cluster=$cluster --save_dir=$SAVED_DATA --ext_annot=$EXT_ANNOT "${FUNCS[@]}"
fi
if [ $? -ne 0 ]; then
  echo "make_report_elements.R has not finished successfully."
  exit 1
fi
Rscript make_report.R --proj=$PROJ --cluster=$cluster --save_dir=$SAVED_DATA --make_single_file=$MAKE_SINGLE_FILE "${FUNCS[@]}"
if [ $? -ne 0 ]; then
  echo "make_report.R has not finished successfully."
  exit 1
fi


for a_func in "${FUNCS[@]}"; do
  case $a_func in
	DEGO_stage)
        Rscript -e rmarkdown::render"('DE-GO-analysis-stagesVS.Rmd',
                 output_file=\"report${PROJ}/data/gv.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	DEGO_1v1)
        Rscript -e rmarkdown::render"('DE-GO-analysis-1v1.Rmd',
                 output_file=\"report${PROJ}/data/1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark_1v1)
        Rscript -e rmarkdown::render"('hallmark-1v1.Rmd',
                 output_file=\"report${PROJ}/data/hallmark_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_1v1)
        Rscript -e rmarkdown::render"('reactome-1v1.Rmd',
                 output_file=\"report${PROJ}/data/reactome_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_1v1)
        Rscript -e rmarkdown::render"('kegg-1v1.Rmd',
                 output_file=\"report${PROJ}/data/kegg_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark_stage)
        Rscript -e rmarkdown::render"('hallmark-stageVS.Rmd',
                 output_file=\"report${PROJ}/data/hallmark_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_stage)
        Rscript -e rmarkdown::render"('reactome-stageVS.Rmd',
                 output_file=\"report${PROJ}/data/reactome_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_stage)
        Rscript -e rmarkdown::render"('kegg-stageVS.Rmd',
                 output_file=\"report${PROJ}/data/kegg_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;


  esac

done
