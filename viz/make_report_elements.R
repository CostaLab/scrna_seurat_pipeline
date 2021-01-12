
# make Rmd report
args <- commandArgs(TRUE)

PROJ = args[1]
CLUSTER = args[2]
SAVE_DIR = args[3]
ANNOTATION_EXTERNAL_FILE_PATH = args[4]
FUNCS = args[-c(1:4)]

# TODO better param processing
if(grepl("--proj",PROJ,fixed=TRUE)) PROJ=gsub("--proj=","",PROJ,fixed=TRUE)
if(grepl("--cluster",CLUSTER,fixed=TRUE)) CLUSTER=gsub("--cluster=","",CLUSTER,fixed=TRUE)
if(grepl("--save_dir",SAVE_DIR,fixed=TRUE)) SAVE_DIR=gsub("--save_dir=","",SAVE_DIR,fixed=TRUE)
if(grepl("--ext_annot",ANNOTATION_EXTERNAL_FILE_PATH,fixed=TRUE)) ANNOTATION_EXTERNAL_FILE_PATH=gsub("--ext_annot=","",ANNOTATION_EXTERNAL_FILE_PATH,fixed=TRUE)


library(Seurat)
library(Hmisc)
library(knitr)
library(kableExtra)
library(dplyr)
library(glue)
library(ggplot2)
library(cowplot)
library(reshape2)
library(urltools)
library(clustree)
library(stringr)
library(digest)
library(openxlsx)
library(ComplexHeatmap)
library(EnhancedVolcano)

# define colour palette
colours = c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99"
)

# define helper functions
`%ni%` <- Negate(`%in%`)

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
    top5.up <- de.data %>% group_by(cluster) %>% top_n(10, avg_logFC) %>%filter(avg_logFC > 0) %>% arrange(-avg_logFC)
    top5.dn <- de.data %>% group_by(cluster) %>% top_n(10, -avg_logFC) %>%filter(avg_logFC < 0) %>% arrange(-avg_logFC)
  } else {
    top5.up <- de.data  %>% top_n(10, avg_logFC) %>%filter(avg_logFC > 0) %>% arrange(-avg_logFC)
    top5.dn <- de.data  %>% top_n(10, -avg_logFC) %>%filter(avg_logFC < 0) %>% arrange(-avg_logFC)
  }
  top.up.dn <- rbind(top5.up, top5.dn)
  top.up.dn$gene <- make.unique(top.up.dn$gene)
  top.up.dn$type = ifelse(top.up.dn$avg_logFC > 0, "positive", "negative")
  top.up.dn$type <- factor(top.up.dn$type, levels = c("positive", "negative"))
  g <- ggplot(data = top.up.dn,
              aes(x = gene, y = avg_logFC, fill = type)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=rev(top.up.dn$gene)) +
    theme_minimal() +
    theme(legend.position="none", axis.text=element_text(size=15)) +
    scale_fill_manual(values = c(positive = "#E41A1C", negative = "#377EB8")) +
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


savedir = SAVE_DIR
cluster = CLUSTER
project = PROJ
funcs = FUNCS
ext_annot_fp = ANNOTATION_EXTERNAL_FILE_PATH

available_clusters = cluster
if(any(grepl("Clusters_harmony",funcs,fixed=TRUE))) available_clusters = c(available_clusters,"harmony_inte_clusters")
if(any(grepl("Clusters_seurat",funcs,fixed=TRUE))) available_clusters = c(available_clusters,"seurat_inte_clusters")
print(available_clusters)
# define/create dirs
report_tables_folder = file.path(paste0('report',PROJ),'tables')
report_plots_folder = file.path(paste0('report',PROJ),'plots')
report_plots_folder_png = file.path(report_plots_folder,'png')
report_plots_folder_pdf = file.path(report_plots_folder,'pdf')

dir.create(report_plots_folder_png, recursive = TRUE)
dir.create(report_plots_folder_pdf, recursive = TRUE)
dir.create(report_tables_folder, recursive = TRUE)

# run necessary generators
if(any(grepl("QC",funcs,fixed=TRUE))) source("1_quality_report_elements.R")
if(any(grepl("DEs",funcs,fixed=TRUE))) source("2_clusters_DEs_elements.R")
if(any(grepl("Clusters",funcs,fixed=TRUE))) source("2_clustering_elements.R")
if(any(grepl("EXT_MARKERS",funcs,fixed=TRUE))) source("3_external_markers_elements.R")
# FIXME possible problem where all term enrichment analysis is on the same place
if(any(grepl("DEGO",funcs,fixed=TRUE))) source("3_DE_GO-analysis_elements.R")
