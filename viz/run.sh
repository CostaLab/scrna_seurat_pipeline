#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' # No Color

PROJ=$1
SAVED_DATA="../save${PROJ}"

FUNCS=(
        --QC
        # --DEs
        # --Clusters
        --Singleton
        # --DEGO
        # --EXT_MARKERS
        # --DEGO_1v1
        # --DEGO_stage
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

cluster="singleton"

#!!!!!!!!!----------------------------------

echo $PROJ
echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p report${PROJ}/data
cp -r ../charts${PROJ}/* report${PROJ}/data
python code_generator.py -c $cluster -cf "conf/config${PROJ}.R" -b "../" -p "$PROJ"
grip --export report${PROJ}/index.md

Rscript make_report.R --proj=$PROJ --cluster=$cluster --save_dir=$SAVED_DATA "${FUNCS[@]}"

# Rscript  -e "rmarkdown::render(
#     'scrna_pipeline_report.Rmd',
#     output_file='report${PROJ}/scrna_report.html',
#     clean=TRUE,
#     params=list(cluster='${cluster}',project='${PROJ}')
#   )" ${SAVED_DATA} "${FUNCS[@]}"

# for a_func in "${FUNCS[@]}"; do
#   case $a_func in
# 	QC)
#         Rscript  -e "rmarkdown::render('1_quality_report.Rmd',
#                  output_file='report${PROJ}/data/data_quality.html',
#                  clean=TRUE)" ${SAVED_DATA}
# 		;;
# 	DEs)
#         Rscript -e "rmarkdown::render('2_clusters_DEs.Rmd',
#                  output_file='report${PROJ}/data/clusters_DEs.html',
#                  clean=TRUE)" ${SAVED_DATA}
# 		;;
# 	Clusters)
#         Rscript -e "rmarkdown::render('2_clustering.Rmd',
#                  output_file='report${PROJ}/data/clusters.html',
#                  clean=TRUE,
#                  params=list(cluster='${cluster}'))" ${SAVED_DATA}
# 		;;
# 	DEGO)
#         Rscript -e "rmarkdown::render('3_DE_GO-analysis.Rmd',
#                  output_file='report${PROJ}/data/dego.html',
#                  clean=TRUE,
#                  params=list(cluster='${cluster}'))" ${SAVED_DATA}
# 		;;
# 	EXT_MARKERS)
#         Rscript -e "rmarkdown::render('3_external_markers.Rmd',
#                  output_file='report${PROJ}/data/external_markers.html',
#                  clean=TRUE,
#                  params=list(cluster='${cluster}'))" ${SAVED_DATA}
# 		;;
#
# 	DEGO_stage)
#         Rscript -e "rmarkdown::render('DE-GO-analysis-stagesVS.Rmd',
#                  output_file='report${PROJ}/data/gv.html',
#                  clean=TRUE,
#                  params=list(cluster='${cluster}'))" ${SAVED_DATA}
# 		;;
# 	DEGO_1v1)
#         Rscript -e "rmarkdown::render('DE-GO-analysis-1v1.Rmd',
#                  output_file='report${PROJ}/data/1vs1.html',
#                  clean=TRUE,
#                  params=list(cluster='${cluster}'))" ${SAVED_DATA}
# 		;;
#   esac
#
# done
