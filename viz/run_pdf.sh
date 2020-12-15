#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' # No Color

FUNCS=(
        QC
        DEs
        Clusters
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
mkdir -p report_pdf/data
cp -r ../charts/* report_pdf/data
python code_generator.py -c $cluster
grip --export report_pdf/index.md



for a_func in "${FUNCS[@]}"; do
  case $a_func in
	QC)
        Rscript -e rmarkdown::render"('1_quality_report.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/data_quality.pdf\",
                 clean=TRUE)"
		;;
	DEs)
        Rscript -e rmarkdown::render"('2_clusters_DEs.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/clusters_DEs.pdf\",
                 clean=TRUE)"
		;;
	Clusters)
        Rscript -e rmarkdown::render"('2_clustering.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/clusters.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	DEGO)
        Rscript -e rmarkdown::render"('3_DE_GO-analysis.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/dego.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	KEGG)
        Rscript -e rmarkdown::render"('3_KEGG.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/KEGG.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark)
        Rscript -e rmarkdown::render"('3_hallmark.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/hallmark.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	Reactome)
        Rscript -e rmarkdown::render"('3_Reactome.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/Reactome.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;



	EXT_MARKERS)
        Rscript -e rmarkdown::render"('3_external_markers.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/external_markers.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	DEGO_stage)
        Rscript -e rmarkdown::render"('DE-GO-analysis-stagesVS.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/gv.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	DEGO_1v1)
        Rscript -e rmarkdown::render"('DE-GO-analysis-1v1.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/1vs1.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark_1v1)
        Rscript -e rmarkdown::render"('hallmark-1v1.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/hallmark_1vs1.pdf\", 
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_1v1)
        Rscript -e rmarkdown::render"('reactome-1v1.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/reactome_1vs1.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_1v1)
        Rscript -e rmarkdown::render"('kegg-1v1.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/kegg_1vs1.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	hallmark_stage)
        Rscript -e rmarkdown::render"('hallmark-stageVS.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/hallmark_stageVS.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;

	reactome_stage)
        Rscript -e rmarkdown::render"('reactome-stageVS.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/reactome_stageVS.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;
	kegg_stage)
        Rscript -e rmarkdown::render"('kegg-stageVS.Rmd', 'pdf_document',
                 output_file=\"report_pdf/data/kegg_stageVS.pdf\",
                 clean=TRUE,
                 params=list(cluster=\"${cluster}\"))"
		;;


  esac


done


