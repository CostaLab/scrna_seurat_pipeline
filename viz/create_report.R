#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(urltools))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(Seurat))
`%ni%` <- Negate(`%in%`)

###===================================FUNCTIONS BEGIN===================================================
run_shell <- function(cmd){
  system(cmd)
}

save_ggplot_formats = function(
  plt, base_plot_dir, plt_name, create_plot_subdir=TRUE,
  formats=c("png","pdf"), units="in", width=20, height=20,
  type="cairo",
  ...
){
  # if theres a plot and basedir
  if(!is.null(base_plot_dir) & !is.null(plt)){
    # for each format
    for(fmt in formats){
      f_path_fmt = file.path(base_plot_dir,paste0(plt_name,".",fmt))
      if(create_plot_subdir) dir.create(file.path(base_plot_dir,fmt),recursive=TRUE,showWarnings=FALSE)
      if(dir.exists(file.path(base_plot_dir,fmt))) f_path_fmt = file.path(base_plot_dir,fmt,paste0(plt_name,".",fmt))
      if(fmt == "png"){
        ggplot2::ggsave(filename=f_path_fmt,plot = plt,device = fmt,units = units,width = width,height = height,type=type,...)
      }else{
        ggplot2::ggsave(filename=f_path_fmt,plot = plt,device = fmt,units = units,width = width,height = height,...)
      }

    }
  }
}

GeneBarPlot <- function(de.data, xlim = NULL, main = NULL) {
  #de.data = cluster.de[[id]]
  #de.data = plot_de
  if (any(colnames(de.data) == "cluster")) {
    top5.up <- de.data %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data %>% group_by(cluster) %>% top_n(10, -avg_log2FC) %>%filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  } else {
    top5.up <- de.data  %>% top_n(10, avg_log2FC) %>%filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data  %>% top_n(10, -avg_log2FC) %>%filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  }
  top.up.dn <- rbind(top5.up, top5.dn)
  top.up.dn$gene <- make.unique(top.up.dn$gene)
  top.up.dn$type = ifelse(top.up.dn$avg_log2FC > 0, "positive", "negative")
  top.up.dn$type <- factor(top.up.dn$type, levels = c("positive", "negative"))
  g <- ggplot(data = top.up.dn,
              aes(x = gene, y = avg_log2FC, fill = type)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=rev(top.up.dn$gene)) +
    theme_minimal() +
    theme(legend.position="none", axis.text=element_text(size=15)) +
    # scale_fill_manual(values = c(positive = "#E41A1C", negative = "#377EB8")) +
    scale_fill_manual(values = c(positive = pos_color, negative = neg_color)) +
    coord_flip()
  if (!is.null(main)) {
    g <- g + ggtitle(main)
  } else {
    g <- g + ggtitle("Average logFC for the top 5 up and top 5 down regulated genes")
  }
  if (!is.null(xlim)) {
    # Coordinates are flipped
    g <- g + ylim(xlim)
  }
  return(g)
}

BarPlot <- function(x, fill, xlab = NULL, ylab = "Cells", legend.title = NULL,
    main = NULL, print.number = TRUE) {
    counts <- as.data.frame(table(x, fill))
    names(counts) <- c("x", "fill", "Cells")
    p <- ggplot(data = counts, aes(x = x, y = Cells, fill = fill)) + geom_bar(stat = "identity",
        position = position_dodge())
    if (print.number) {
        p <- p + geom_text(aes(label = Cells), vjust = 1.6, color = "white",
            position = position_dodge(0.9), size = 3.5)
    }
    if (!is.null(xlab)) {
        p <- p + xlab(xlab)
    }
    if (!is.null(ylab)) {
        p <- p + ylab(ylab)
    }
    if (!is.null(legend.title)) {
        p <- p + labs(fill = legend.title)
    }
    if (!is.null(main)) {
        p <- p + ggtitle(main)
    }
    return(p)
}
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}

###===================================FUNCTIONS END===================================================

###===================================PARAMETERS BEGIN===================================================
AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-a", "--author"), type="character", default="Costalab",
                       help="Author display in report [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-m", "--make_element"), type="character", default="FALSE",
                       help="make_element [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-s", "--savedir"), type="character", default="save",
                       help="savedir [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-c", "--configfile"), type="character", default="conf/config.R",
                       help="configfile [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-o", "--report_dir"), type="character", default="report",
                       help="report_dir [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-r", "--charts_dir"), type="character", default="charts",
                       help="charts_dir [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-p", "--project"), type="character", default="",
                       help="charts_dir [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-e", "--externalfile"), type="character", default="./external/Human_and_mouse_cell_markers-Markers.tsv",
                       help="report_dir [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-d", "--defaultclsuters"), type="character", default="seurat_clusters",
                       help="cluster slots: \nremoved_clusters\nremove_reclusters\nmerged_clusters\nannotation\nsingleton [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-j", "--planOfreport"), type="character", default="c()",
                       help="plan of report [default %default]",
                       metavar="character")
  parser <- add_option(parser, c("-l", "--singlefile"), type="character", default="FALSE",
                       help="generate single html file [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("-i", "--indexonly"), type="character", default="FALSE",
                       help="only generate index.html [default %default]",
                       metavar="character")



  return(parser)
}
parser <- AllOptions()
#debug(parse_args)
pa <- parse_args(parser)
#print(pa)

PROJECT           = pa$project
AUTHOR            = pa$author
MAKE_ELEMENT      = pa$make_element
SAVE_DIR          = pa$savedir
CONFIGFILE        = pa$configfile
REPORTDIR         = pa$report_dir
CHARTSDIR         = pa$charts_dir
EXTERNALFILE      = pa$externalfile
DEFAULTCLUSTERS   = pa$defaultclsuters
GEN_SINGLE_FILE   = pa$singlefile
EXEC_PLAN         = jsonlite::fromJSON(pa$planOfreport)
INDEX_ONLY        = pa$indexonly

#REPORTDIR     = paste0(REPORTDIR, "_", PROJECT)
#print(REPORTDIR)
dir.create(file.path(REPORTDIR, "data"), recursive=T)


py_exe_list = paste0(EXEC_PLAN, collapse="','")
py_exe_list = paste0("['", py_exe_list, "']")
###===================================PARAMETERS END===================================================

message("Project:", PROJECT)
cat(paste("Use cluster slot", red(DEFAULTCLUSTERS), "\n"))
##1. copy excels to report/data
message("copying excel files from ", CHARTSDIR, "\n")
file.copy(file.path(CHARTSDIR, list.files(CHARTSDIR)), glue("{REPORTDIR}/data"), overwrite = T, copy.date=T)
dir.create(glue("{REPORTDIR}/../viz"))
file.copy(file.path(getwd(), "viz", list.files("viz")), glue("{REPORTDIR}/../viz"), overwrite = T, recursive=T, copy.date=T)

code_gene_file <- file.path(REPORTDIR, "..", "viz/code_generator.py")

viz_path <- file.path(REPORTDIR, "..", "viz")

code_generate_cmd=glue("python {code_gene_file} ",
                       "-c {DEFAULTCLUSTERS} ",
                       "-cf '{CONFIGFILE}' ",
                       "--save_dir '{SAVE_DIR}' ",
                       "--output_dir '{REPORTDIR}' ",
                       "--proj_tag {PROJECT} ",
                       "--executing_list \"{py_exe_list}\"")


##2. code generate
run_shell(code_generate_cmd)
##3. generate index.html from template
run_shell(glue("grip --export {REPORTDIR}/index.md"))

#stop(INDEX_ONLY)
if(INDEX_ONLY == "TRUE"){
  cat(red("======Only generate index.html=====\n"))
  quit(save='no')
}

report_data_folder  = file.path(REPORTDIR, "data")
report_tables_folder = file.path(REPORTDIR,'tables')
report_plots_folder = file.path(REPORTDIR,'plots')
report_plots_folder_png = file.path(report_plots_folder,'png')
report_plots_folder_pdf = file.path(report_plots_folder,'pdf')
savedir = SAVE_DIR
cluster = DEFAULTCLUSTERS
project = PROJECT
funcs = EXEC_PLAN
ext_annot_fp = EXTERNALFILE



##4. Make_report element
  if (MAKE_ELEMENT == "TRUE"){
  source(CONFIGFILE)
  cluster_viridis_opt = ifelse(
    any(grepl("cluster_color_option",names(viz_conf),fixed = TRUE)),
    viz_conf[["cluster_color_option"]], # Config option
    "D" # Default
  )

  # replicate_colors
  replicates_viridis_opt = ifelse(
    any(grepl("replicate_color_option",names(viz_conf),fixed = TRUE)),
    viz_conf[["replicate_color_option"]],
    "C"
  )
  neg_color = ifelse(
    any(grepl("neg_color",names(viz_conf),fixed = TRUE)),
    viz_conf[["neg_color"]],
    '#51C3CC' # colorBlindness::Blue2DarkOrange12Steps[2]
  )

  # divergent color for positive values
  pos_color = ifelse(
    any(grepl("pos_color",names(viz_conf),fixed = TRUE)),
    viz_conf[["pos_color"]],
    '#CC5800' #rev(colorBlindness::Blue2DarkOrange12Steps)[2]
  )

  # divergent color for base values
  base_color = ifelse(
    any(grepl("base_color",names(viz_conf),fixed = TRUE)),
    viz_conf[["base_color"]],
    "lightgrey"
  )

  # divergent palette for neg to pos values
  neg_pos_divergent_palette = ifelse(
    any(grepl("neg_pos_divergent_palette",names(viz_conf),fixed = TRUE)),
    viz_conf[["neg_pos_divergent_palette"]],
    # colorBlindness::Blue2DarkOrange12Steps
    c('#1E8E99','#51C3CC','#99F9FF','#B2FCFF','#CCFEFF','#E5FFFF','#FFE5CC','#FFCA99','#FFAD65','#FF8E32','#CC5800','#993F00')
  )

  zero_pos_divergent_colors = c(base_color,pos_color)

  override_colours = c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
    "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99"
  )



  available_clusters = DEFAULTCLUSTERS
  if("Clusters_harmony" %in% EXEC_PLAN){
    available_clusters <- c(available_clusters, "harmony_inte_clusters")
  }
  if("Clusters_seurat" %in% EXEC_PLAN){
    available_clusters = c(available_clusters,"seurat_inte_clusters")
  }


  #paste("report_data_folder", blue(report_data_folder))
  #paste("report_data_folder"blue(report_table_folder))

  dir.create(report_plots_folder_png, recursive = TRUE)
  dir.create(report_plots_folder_pdf, recursive = TRUE)
  dir.create(report_tables_folder, recursive = TRUE)

  # run necessary generators
  if("QC" %in% EXEC_PLAN) {
    cat(paste(date(), green(" Element: "), red("QC"), "\n"))
    source(glue("{viz_path}/1_quality_report_elements.R"))
  }
  if(any(grepl("^AmbientRNA",funcs,fixed=TRUE))){
    cat(paste(date(), green(" Element: "), red("AmbientRNA"), "\n"))
    source(glue("{viz_path}/ambientRNA_viz_elements.R"))
  }
  if(any(grepl("^soupxAmbientRNA",funcs,fixed=TRUE))){
    cat(paste(date(), green(" Element: "), red("soupxAmbientRNA"), "\n"))
    source(glue("{viz_path}/ambientRNA_soupx_viz_elements.R"))
  }
  if("DEs"%in% EXEC_PLAN) {
    cat(paste(date(), green(" Element: "), red("DEs"), "\n"))
    source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  }
  if(any(grepl("Clusters_", EXEC_PLAN, fixed=T))){
    cat(paste(date(), green(" Element: "), red("Clusters_?"), "\n"))
    source(glue("{viz_path}/2_batch_clustering_elements.R"))
    source(glue("{viz_path}/2_clustering_elements.R"))
  }
  if("Clusters" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("Clusters"), "\n"))
    source(glue("{viz_path}/2_clustering_elements.R"))
  }
  if("Singleton" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("Singleton"), "\n"))
    source(glue("{viz_path}/2_clustering_elements.R"))
    source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  }
  if("EXT_MARKERS" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("EXT_MARKERS"), "\n"))
    source(glue("{viz_path}/3_external_markers_elements.R"))
  }
  # FIXME possible problem where all term enrichment analysis is on the same place
  if("DEGO" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("DEGO"), "\n"))
    source(glue("{viz_path}/3_DE_GO-analysis_elements.R"))
  }

}

##5. Produce Report
render_func = function(rmd_input_filename, output_filename){
  rmarkdown::render(
    rmd_input_filename,
    output_file=file.path(REPORTDIR,"data",output_filename),
    output_format=c("html_document"),
    clean=TRUE,
    params=list(
      cluster=DEFAULTCLUSTERS,
      project=PROJECT,
      savedir=SAVE_DIR,
      funcs=EXEC_PLAN,
      report_data_folder=report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author=AUTHOR
    )
  )
}

for(i in EXEC_PLAN){
  cat(paste(date(), blue(" Generating: "), red(i), "\n"))

  if(grepl("QC",i,fixed=TRUE)) render_func(glue("{viz_path}/1_quality_report.Rmd"),"data_quality")
  if(grepl("AmbientRNA",i,fixed=TRUE)) render_func(glue("{viz_path}/ambientRNA_viz.Rmd"),"ambient_rna")
  if(grepl("^soupxAmbientRNA",i,fixed=TRUE)) render_func(glue("{viz_path}/ambientRNA_soupx_viz.Rmd"),"ambient_rna_soupx")
  if(grepl("DEs",i,fixed=TRUE)) render_func(glue("{viz_path}/2_clusters_DEs.Rmd"),"clusters_DEs")
  if(grepl("Clusters",i,fixed=TRUE) | grepl("Singleton",i,fixed=TRUE)) render_func(glue("{viz_path}/2_clustering.Rmd"),"clusters")
  if(grepl("Clusters_harmony",i,fixed=TRUE)) render_func(glue("{viz_path}/2_clustering_harmony.Rmd"),"clusters_harmony")
  if(grepl("Clusters_seurat",i,fixed=TRUE)) render_func(glue("{viz_path}/2_clustering_seurat.Rmd"),"clusters_seurat")
  if(grepl("EXT_MARKERS",i,fixed=TRUE)) render_func(glue("{viz_path}/3_external_markers.Rmd"),"external_markers")
  if(grepl("DEGO",i,fixed=TRUE)) render_func(glue("{viz_path}/3_DE_GO-analysis.Rmd"),"dego")
  if(grepl("KEGG",i,fixed=TRUE)) render_func(glue("{viz_path}/3_KEGG.Rmd"),"KEGG")
  if(grepl("progeny",i,fixed=TRUE)) render_func(glue("{viz_path}/3_progeny.Rmd"),"progeny")
  if(grepl("hallmark",i,fixed=TRUE)) render_func(glue("{viz_path}/3_hallmark.Rmd"),"hallmark")
  if(grepl("Reactome",i,fixed=TRUE)) render_func(glue("{viz_path}/3_Reactome.Rmd"),"Reactome")
  if(grepl("DEGO_stage",i,fixed=TRUE)) render_func(glue("{viz_path}/DE-GO-analysis-stagesVS.Rmd"),"gv")
  if(grepl("DEGO_1v1",i,fixed=TRUE)) render_func(glue("{viz_path}/DE-GO-analysis-1v1.Rmd"),"1vs1")
  if(grepl("hallmark_1v1",i,fixed=TRUE)) render_func(glue("{viz_path}/hallmark-1v1.Rmd"),"hallmark_1vs1")
  if(grepl("reactome_1v1",i,fixed=TRUE)) render_func(glue("{viz_path}/reactome-1v1.Rmd"),"reactome_1vs1")
  if(grepl("kegg_1v1",i,fixed=TRUE)) render_func(glue("{viz_path}/kegg-1v1.Rmd"),"kegg_1vs1")
  if(grepl("hallmark_stage",i,fixed=TRUE)) render_func(glue("{viz_path}/hallmark-stageVS.Rmd"),"hallmark_stageVS")
  if(grepl("reactome_stage",i,fixed=TRUE)) render_func(glue("{viz_path}/reactome-stageVS.Rmd"),"reactome_stageVS")
  if(grepl("progeny_stage",i,fixed=TRUE)) render_func(glue("{viz_path}/progeny-stageVS.Rmd"),"progeny_stageVS")
  if(grepl("kegg_stage",i,fixed=TRUE)) render_func(glue("{viz_path}/kegg-stageVS.Rmd"),"kegg_stageVS")
  if(grepl("intUMAPs",i,fixed=TRUE)) render_func(glue("{viz_path}/interactive_UMAPs.Rmd"),"interactive_UMAPs")
}

if(GEN_SINGLE_FILE){
  # generate report including everything in a single file
  rmarkdown::render(
    glue('{viz_path}/scrna_pipeline_report.Rmd'),
    output_file=file.path(REPORTDIR,'scrna_report'),
    output_format=c("html_document"),
    clean=TRUE,
    params=list(
      cluster=DEFAULTCLUSTERS,
      project=PROJECT,
      savedir=SAVE_DIR,
      funcs=EXEC_PLAN,
      report_data_folder=report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author=AUTHOR
    )
  )
}

