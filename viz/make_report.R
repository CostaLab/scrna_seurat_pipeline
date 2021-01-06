
# make Rmd report
args <- commandArgs(TRUE)

PROJ = args[1]
CLUSTER = args[2]
SAVE_DIR = args[3]
FUNCS = args[-c(1:3)]

# TODO better param processing
if(grepl("--proj",PROJ,fixed=TRUE)) PROJ=gsub("--proj=","",PROJ,fixed=TRUE)
if(grepl("--cluster",CLUSTER,fixed=TRUE)) CLUSTER=gsub("--cluster=","",CLUSTER,fixed=TRUE)
if(grepl("--save_dir",SAVE_DIR,fixed=TRUE)) SAVE_DIR=gsub("--save_dir=","",SAVE_DIR,fixed=TRUE)

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



# report_data = as.list(list.files(file.path(paste0('report',PROJ),'data')))
report_data_folder = file.path(paste0('report',PROJ),'data')

rmarkdown::render(
  'scrna_pipeline_report.Rmd',
  output_file=file.path(paste0('report',PROJ),'scrna_report'),
  output_format=c("html_document"),
  clean=TRUE,
  params=list(
    cluster=CLUSTER,
    project=PROJ,
    savedir=SAVE_DIR,
    funcs=FUNCS,
    report_data_folder=report_data_folder,
    report_tables_folder = file.path(paste0('report',PROJ),'tables'),
    report_plots_folder  = file.path(paste0('report',PROJ),'plots'),
    report_plots_folder_png = file.path(paste0('report',PROJ),'plots','png'),
    report_plots_folder_pdf = file.path(paste0('report',PROJ),'plots','pdf')
  )
)

render_func = function(rmd_input_filename, output_filename){
  rmarkdown::render(
    rmd_input_filename,
    output_file=file.path(paste0('report',PROJ),"data",output_filename),
    output_format=c("html_document"),
    clean=TRUE,
    params=list(
      cluster=CLUSTER,
      project=PROJ,
      savedir=SAVE_DIR,
      funcs=FUNCS,
      report_data_folder=report_data_folder,
      report_tables_folder = file.path(paste0('report',PROJ),'tables'),
      report_plots_folder  = file.path(paste0('report',PROJ),'plots'),
      report_plots_folder_png = file.path(paste0('report',PROJ),'plots','png'),
      report_plots_folder_pdf = file.path(paste0('report',PROJ),'plots','pdf')
    )
  )
}

for(i in FUNCS){
  if(grepl("QC",i,fixed=TRUE)) render_func("1_quality_report.Rmd","data_quality")
  if(grepl("DEs",i,fixed=TRUE)) render_func("2_clusters_DEs.Rmd","clusters_DEs")
  if(grepl("Clusters",i,fixed=TRUE)) render_func("2_clustering.Rmd","clusters")
  if(grepl("Clusters_harmony",i,fixed=TRUE)) render_func("2_clustering_harmony.Rmd","clusters_harmony")
  if(grepl("Clusters_seurat",i,fixed=TRUE)) render_func("2_clustering_seurat.Rmd","clusters_seurat")
  if(grepl("DEGO",i,fixed=TRUE)) render_func("3_DE_GO-analysis.Rmd","dego")
  if(grepl("KEGG",i,fixed=TRUE)) render_func("3_KEGG.Rmd","KEGG")
  if(grepl("hallmark",i,fixed=TRUE)) render_func("3_hallmark.Rmd","hallmark")
  if(grepl("Reactome",i,fixed=TRUE)) render_func("3_Reactome.Rmd","Reactome")
  if(grepl("EXT_MARKERS",i,fixed=TRUE)) render_func("3_external_markers.Rmd","external_markers")
  if(grepl("DEGO_stage",i,fixed=TRUE)) render_func("DE-GO-analysis-stagesVS.Rmd","gv")
  if(grepl("DEGO_1v1",i,fixed=TRUE)) render_func("DE-GO-analysis-1v1.Rmd","1vs1")
  if(grepl("hallmark_1v1",i,fixed=TRUE)) render_func("hallmark-1v1.Rmd","hallmark_1vs1")
  if(grepl("reactome_1v1",i,fixed=TRUE)) render_func("reactome-1v1.Rmd","reactome_1vs1")
  if(grepl("kegg_1v1",i,fixed=TRUE)) render_func("kegg-1v1.Rmd","kegg_1vs1")
  if(grepl("hallmark_stage",i,fixed=TRUE)) render_func("hallmark-stageVS.Rmd","hallmark_stageVS")
  if(grepl("reactome_stage",i,fixed=TRUE)) render_func("reactome-stageVS.Rmd","reactome_stageVS")
  if(grepl("kegg_stage",i,fixed=TRUE)) render_func("kegg-stageVS.Rmd","kegg_stageVS")
}
