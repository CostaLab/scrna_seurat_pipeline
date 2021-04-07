
# make Rmd report
args <- commandArgs(TRUE)

PROJ = args[1]
CLUSTER = args[2]
SAVE_DIR = args[3]
GEN_SINGLE_FILE = args[4]
AUTHOR = args[5]
OUTPUT_DIR = args[6]
FUNCS = args[-c(1:6)]

library(Seurat)
library(celda)
library(ggplot2)

# TODO better param processing
if(grepl("--proj",PROJ,fixed=TRUE)) PROJ=gsub("--proj=","",PROJ,fixed=TRUE)
if(grepl("--cluster",CLUSTER,fixed=TRUE)) CLUSTER=gsub("--cluster=","",CLUSTER,fixed=TRUE)
if(grepl("--save_dir",SAVE_DIR,fixed=TRUE)) SAVE_DIR=gsub("--save_dir=","",SAVE_DIR,fixed=TRUE)
if(grepl("--make_single_file",GEN_SINGLE_FILE,fixed=TRUE)){
  GEN_SINGLE_FILE=as.logical(gsub("--make_single_file=","",GEN_SINGLE_FILE,fixed=TRUE))
}else{
  GEN_SINGLE_FILE=FALSE
}
if(grepl("--author",AUTHOR,fixed=TRUE)) AUTHOR=gsub("--author=","",AUTHOR,fixed=TRUE)
if(grepl("--output_dir",OUTPUT_DIR,fixed=TRUE)) OUTPUT_DIR=gsub("--output_dir=","",OUTPUT_DIR,fixed=TRUE)

# TODO
# might need system packages to create pdf version
# texlive-latex-base
# texlive-latex-recommended
# texlive-fonts-recommended
# texlive-latex-extra

# helper function
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}


report_data_folder = file.path(OUTPUT_DIR,'data')
report_tables_folder = file.path(OUTPUT_DIR,'tables')
report_plots_folder = file.path(OUTPUT_DIR,'plots')
report_plots_folder_png = file.path(OUTPUT_DIR,'plots','png')
report_plots_folder_pdf = file.path(OUTPUT_DIR,'plots','pdf')

render_func = function(rmd_input_filename, output_filename){
  rmarkdown::render(
    rmd_input_filename,
    output_file=file.path(OUTPUT_DIR,"data",output_filename),
    output_format=c("html_document"),
    clean=TRUE,
    params=list(
      cluster=CLUSTER,
      project=PROJ,
      savedir=SAVE_DIR,
      funcs=FUNCS,
      report_data_folder=report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author=AUTHOR
    )
  )
}

for(i in FUNCS){
<<<<<<< HEAD
  if(grepl("QC",i,fixed=TRUE)) render_func("1_quality_report.Rmd","data_quality")
  if(grepl("AmbientRNA",i,fixed=TRUE)) render_func("ambientRNA_viz.Rmd","ambient_rna")
  if(grepl("DEs",i,fixed=TRUE)) render_func("2_clusters_DEs.Rmd","clusters_DEs")
  if(grepl("Clusters",i,fixed=TRUE) | grepl("Singleton",i,fixed=TRUE)) render_func("2_clustering.Rmd","clusters")
  if(grepl("Clusters_harmony",i,fixed=TRUE)) render_func("2_clustering_harmony.Rmd","clusters_harmony")
  if(grepl("Clusters_seurat",i,fixed=TRUE)) render_func("2_clustering_seurat.Rmd","clusters_seurat")
  if(grepl("EXT_MARKERS",i,fixed=TRUE)) render_func("3_external_markers.Rmd","external_markers")
  if(grepl("DEGO",i,fixed=TRUE)) render_func("3_DE_GO-analysis.Rmd","dego")
  if(grepl("KEGG",i,fixed=TRUE)) render_func("3_KEGG.Rmd","KEGG")
  if(grepl("hallmark",i,fixed=TRUE)) render_func("3_hallmark.Rmd","hallmark")
  if(grepl("Reactome",i,fixed=TRUE)) render_func("3_Reactome.Rmd","Reactome")
  if(grepl("DEGO_stage",i,fixed=TRUE)) render_func("DE-GO-analysis-stagesVS.Rmd","gv")
  if(grepl("DEGO_1v1",i,fixed=TRUE)) render_func("DE-GO-analysis-1v1.Rmd","1vs1")
  if(grepl("hallmark_1v1",i,fixed=TRUE)) render_func("hallmark-1v1.Rmd","hallmark_1vs1")
  if(grepl("reactome_1v1",i,fixed=TRUE)) render_func("reactome-1v1.Rmd","reactome_1vs1")
  if(grepl("kegg_1v1",i,fixed=TRUE)) render_func("kegg-1v1.Rmd","kegg_1vs1")
  if(grepl("hallmark_stage",i,fixed=TRUE)) render_func("hallmark-stageVS.Rmd","hallmark_stageVS")
  if(grepl("reactome_stage",i,fixed=TRUE)) render_func("reactome-stageVS.Rmd","reactome_stageVS")
  if(grepl("kegg_stage",i,fixed=TRUE)) render_func("kegg-stageVS.Rmd","kegg_stageVS")
  if(grepl("intUMAPs",i,fixed=TRUE)) render_func("interactive_UMAPs.Rmd","interactive_UMAPs")
=======
  if(grepl("QC",i,fixed=TRUE)) render_func("viz/1_quality_report.Rmd","data_quality")
  if(grepl("DEs",i,fixed=TRUE)) render_func("viz/2_clusters_DEs.Rmd","clusters_DEs")
  if(grepl("Clusters",i,fixed=TRUE) | grepl("Singleton",i,fixed=TRUE)) render_func("viz/2_clustering.Rmd","clusters")
  if(grepl("Clusters_harmony",i,fixed=TRUE)) render_func("viz/2_clustering_harmony.Rmd","clusters_harmony")
  if(grepl("Clusters_seurat",i,fixed=TRUE)) render_func("viz/2_clustering_seurat.Rmd","clusters_seurat")
  if(grepl("EXT_MARKERS",i,fixed=TRUE)) render_func("viz/3_external_markers.Rmd","external_markers")
  if(grepl("DEGO",i,fixed=TRUE)) render_func("viz/3_DE_GO-analysis.Rmd","dego")
  if(grepl("KEGG",i,fixed=TRUE)) render_func("viz/3_KEGG.Rmd","KEGG")
  if(grepl("progeny",i,fixed=TRUE)) render_func("viz/3_progeny.Rmd","progeny")
  if(grepl("hallmark",i,fixed=TRUE)) render_func("viz/3_hallmark.Rmd","hallmark")
  if(grepl("Reactome",i,fixed=TRUE)) render_func("viz/3_Reactome.Rmd","Reactome")
  if(grepl("DEGO_stage",i,fixed=TRUE)) render_func("viz/DE-GO-analysis-stagesVS.Rmd","gv")
  if(grepl("DEGO_1v1",i,fixed=TRUE)) render_func("viz/DE-GO-analysis-1v1.Rmd","1vs1")
  if(grepl("hallmark_1v1",i,fixed=TRUE)) render_func("viz/hallmark-1v1.Rmd","hallmark_1vs1")
  if(grepl("reactome_1v1",i,fixed=TRUE)) render_func("viz/reactome-1v1.Rmd","reactome_1vs1")
  if(grepl("kegg_1v1",i,fixed=TRUE)) render_func("viz/kegg-1v1.Rmd","kegg_1vs1")
  if(grepl("hallmark_stage",i,fixed=TRUE)) render_func("viz/hallmark-stageVS.Rmd","hallmark_stageVS")
  if(grepl("reactome_stage",i,fixed=TRUE)) render_func("viz/reactome-stageVS.Rmd","reactome_stageVS")
  if(grepl("progeny_stage",i,fixed=TRUE)) render_func("viz/progeny-stageVS.Rmd","progeny_stageVS")
  if(grepl("kegg_stage",i,fixed=TRUE)) render_func("viz/kegg-stageVS.Rmd","kegg_stageVS")
  if(grepl("intUMAPs",i,fixed=TRUE)) render_func("viz/interactive_UMAPs.Rmd","interactive_UMAPs")
>>>>>>> master
}


if(GEN_SINGLE_FILE){
  # generate report including everything in a single file
  rmarkdown::render(
    'viz/scrna_pipeline_report.Rmd',
    output_file=file.path(OUTPUT_DIR,'scrna_report'),
    output_format=c("html_document"),
    clean=TRUE,
    params=list(
      cluster=CLUSTER,
      project=PROJ,
      savedir=SAVE_DIR,
      funcs=FUNCS,
      report_data_folder=report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author=AUTHOR
    )
  )
}
