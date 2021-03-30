#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' # No Color

PROJ=$1
SAVED_DATA="../save${PROJ}"

FUNCS=(
        #--QC
        --Ambient_RNA
        #--Singleton
        #--DEs
        # QC
        #Clusters
        #Clusters_harmony
        #Clusters_seurat
        DEGO
        KEGG
        hallmark
        Reactome
        EXT_MARKERS
        DEGO_stage
        hallmark_stage
        reactome_stage
        kegg_stage
        DEGO_1v1
        hallmark_1v1
        reactome_1v1
        kegg_1v1
	#GET_DATA
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

cluster="annotation"

#!!!!!!!!!----------------------------------

echo $PROJ
echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p report/data
cp -r -p ../charts/* report/data
python3 code_generator.py -c $cluster
grip --export report/index.md



for a_func in "${FUNCS[@]}"; do
  case $a_func in
	QC)
        Rscript -e rmarkdown::render"('1_quality_report.Rmd',
                 output_file=\"report/data/data_quality.html\",
                 clean=TRUE)"
	        ;;
	Ambient_RNA)
	Rscript -e rmarkdown::render"('ambientRNA_viz.Rmd',
                 output_file=\"report/data/ambient_rna.html\",
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

	Clusters_harmony)
        Rscript -e rmarkdown::render"('2_clustering_harmony.Rmd',
                 output_file=\"report/data/clusters_harmony.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	Clusters_seurat)
        Rscript -e rmarkdown::render"('2_clustering_seurat.Rmd',
                 output_file=\"report/data/clusters_seurat.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	DEGO)
        Rscript -e rmarkdown::render"('3_DE_GO-analysis.Rmd',
                 output_file=\"report/data/dego.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	KEGG)
        Rscript -e rmarkdown::render"('3_KEGG.Rmd',
                 output_file=\"report/data/KEGG.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark)
        Rscript -e rmarkdown::render"('3_hallmark.Rmd',
                 output_file=\"report/data/hallmark.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	Reactome)
        Rscript -e rmarkdown::render"('3_Reactome.Rmd',
                 output_file=\"report/data/Reactome.html\",
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

	hallmark_1v1)
        Rscript -e rmarkdown::render"('hallmark-1v1.Rmd',
                 output_file=\"report/data/hallmark_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_1v1)
        Rscript -e rmarkdown::render"('reactome-1v1.Rmd',
                 output_file=\"report/data/reactome_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_1v1)
        Rscript -e rmarkdown::render"('kegg-1v1.Rmd',
                 output_file=\"report/data/kegg_1vs1.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark_stage)
        Rscript -e rmarkdown::render"('hallmark-stageVS.Rmd',
                 output_file=\"report/data/hallmark_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_stage)
        Rscript -e rmarkdown::render"('reactome-stageVS.Rmd',
                 output_file=\"report/data/reactome_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_stage)
        Rscript -e rmarkdown::render"('kegg-stageVS.Rmd',
                 output_file=\"report/data/kegg_stageVS.html\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;


  esac


done


Rscript make_report.R --proj=$PROJ --cluster=$cluster --save_dir=$SAVED_DATA "${FUNCS[@]}"
