#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' # No Color

FUNCS=(
        QC 
        DEs
        Clusters
        DEGO
        EXT_MARKERS
        DEGO_1v1
        DEGO_stage 
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

cluster="seurat_clusters"

#!!!!!!!!!----------------------------------


echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p report/data
cp -r ../charts/* report/data
python code_generator.py -c $cluster
grip --export report/index.md



for a_func in "${FUNCS[@]}"; do
  case $a_func in
	QC)
        Rscript -e rmarkdown::render"('1_quality_report.Rmd',
                 output_file=\"report/data/data_quality.html\",
                 clean=TRUE)"
		;;
	DEs)
        Rscript -e rmarkdown::render"('2_clusters_DEs.Rmd',
                 output_file=\"report/data/clusters_DEs.html\",
                 clean=TRUE)"
		;;
	Clusters)
        Rscript -e rmarkdown::render"('2_clustering.Rmd',
                 output_file=\"report/data/clusters.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	DEGO)
        Rscript -e rmarkdown::render"('3_DE_GO-analysis.Rmd',
                 output_file=\"report/data/dego.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	EXT_MARKERS)
        Rscript -e rmarkdown::render"('3_external_markers.Rmd',
                 output_file=\"report/data/external_markers.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	DEGO_stage)
        Rscript -e rmarkdown::render"('DE-GO-analysis-stagesVS.Rmd',
                 output_file=\"report/data/gv.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	DEGO_1v1)
        Rscript -e rmarkdown::render"('DE-GO-analysis-1v1.Rmd',
                 output_file=\"report/data/1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
  esac

done


