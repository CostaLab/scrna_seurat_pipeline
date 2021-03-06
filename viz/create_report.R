#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))      ## Options
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
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DOSE))


`%ni%` <- Negate(`%in%`)

###===================================FUNCTIONS BEGIN===================================================
run_shell <- function(cmd){
  system(cmd)
}

save_ggplot_formats = function(
  plt, base_plot_dir, plt_name, create_plot_subdir=TRUE,
  formats=c("png","pdf"), units="in", width=20, height=20,
  type="cairo",
  plot_obj="ggplot",
  ...
){
  # if theres a plot and basedir
  if(!is.null(base_plot_dir) & !is.null(plt)){
    # for each format
    for(fmt in formats){
      f_path_fmt = file.path(base_plot_dir,paste0(plt_name,".",fmt))
      if(create_plot_subdir) dir.create(file.path(base_plot_dir,fmt),recursive=TRUE,showWarnings=FALSE)
      if(dir.exists(file.path(base_plot_dir,fmt))) f_path_fmt = file.path(base_plot_dir,fmt,paste0(plt_name,".",fmt))
      if(plot_obj == "ggplot"){
        if(fmt == "png"){
          ggplot2::ggsave(filename=f_path_fmt,plot = plt,device = fmt,units = units,width = width,height = height,type=type,...)
        }else{
          ggplot2::ggsave(filename=f_path_fmt,plot = plt,device = fmt,units = units,width = width,height = height,...)
        }
      }else
        if(fmt == "png"){
          png(filename=f_path_fmt,units = units,width = width,height = height, type=type)
          draw(plt, ...)
          dev.off()
        }else{
          pdf(file=f_path_fmt, width = width,height = height)
          draw(plt, ...)
          dev.off()
        }
    }
  }
}

build_cluster_info <- function(scrna){
  if(is.null(scrna)){
    return(DEFAULTCLUSTERS)
  }
  #INTEGRATION_OPTION
  if(DEFAULTCLUSTERS == "seurat_clusters"){
    len <- length(scrna@tools$parameter)
    resol <- scrna@tools$parameter[[len]][["clusterresolution"]]
    #message(resol)
    #INTEGRATION_OPTION
    #message(DEFAULTCLUSTERS)
    #message(INTEGRATION_OPTION)
    st <- glue("{DEFAULTCLUSTERS}, resolution: {resol}, integration option: {INTEGRATION_OPTION}")
    return(st)
  }else{
    st <- glue("{DEFAULTCLUSTERS}, integration option: {INTEGRATION_OPTION}")
    return(st)
  }
  return(DEFAULTCLUSTERS)
}


GeneBarPlot <- function(de.data, xlim = NULL, main = NULL) {
  #de.data = cluster.de[[id]]
  #de.data = plot_de
  if("avg_logFC" %in% names(de.data)){ ## compatible for seurat3
    de.data$avg_log2FC <- de.data$avg_logFC/log(2)
  }
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

is_contigous_true_df <- function(is_sigs){
  ret_df <- data.frame(keep=FALSE, avgIdx=-1)
  if(any(is_sigs) & table(is_sigs)['TRUE'] == 1){
        ret_df$keep=TRUE
        ret_df$avgIdx = which(is_sigs == TRUE)
        return(ret_df)
  }
  return(ret_df)
}

pw_bar_plots <- function(pw_list){
  term_plot_list = lapply(
                          pw_list,
                          function(x){
                            df = x@result
                            if (is.list(df) && length(df)==0){
                              log_m = as.data.frame(list())
                              return(log_m)
                            }
                            log_m = as.data.frame(-log10(df$p.adjust))
                            log_m$names = as.factor(sapply(df$Description, function(y){
                                                             y <- as.character(trimws(y))
                                                             return(y) }))
                            log_m$show_names = as.factor(sapply(df$Description, function(y){
                                                                  y <- as.character(trimws(y))
                                                                  y <- ifelse(nchar(y)<=33,  y, paste0(substr(y, 1, 30), "..."))
                                                                  return(y) }))
                            log_m <- log_m[order(log_m[,1],decreasing = TRUE),]
                            showCatetermry = min(length(log_m[,1]), 10)
                            log_m <- log_m[1:showCatetermry, ]
                            log_m <- log_m[order(log_m[,1],decreasing = FALSE),]
                            return(log_m)
                          }
  )

  ### direction genes plot
  plots <- lapply(
                  seq_along(term_plot_list),
                  function(y, i) {
                    col <- y[[i]]
                    if(length(col) == 0)
                      return(NULL)
                    ggplot(col, aes(reorder(x=col[,2], col[,1]), y=col[,1])) +
                      #ggplot(col, aes(reorder(x=col[,2], col[,1]), y=col[,1])) +
                      geom_bar(stat="identity", fill="#3399CC", color="grey50") +
                      ggtitle(paste(names(y)[i])) +
                      theme_minimal() +
                      theme(axis.text.y  = element_text(size=20), axis.title.y=element_blank(),axis.ticks.y=element_blank()) +
                      labs(y = "") +
                      scale_y_continuous(name="-log10(p-value)") +
                      scale_x_discrete(breaks = col[,2], labels = col[,3]) +
                      coord_flip()
                  },
                  y=term_plot_list
  )
  plots <-Filter(Negate(is.null), plots)
  return(plots)
}



pw_mtx_create <- function(df_list){
  union_TermID <- Reduce(union, lapply(df_list, function(x) x$ID))

  union_df <- do.call(rbind, df_list)
  union_df$GeneRatio <- 0
  union_df$BgRatio <- 0
  union_df$pvalue <- 1
  union_df$p.adjust <- 1
  union_df$qvalue <- 1
  union_df$geneID <- ""
  union_df$Count <- 0

  union_df <- union_df[!duplicated(union_df$ID), ]
  rownames(union_df) <- union_df$ID
  for(x in 1:length(df_list)){
    rest_ids <- setdiff(rownames(union_df), rownames(df_list[[x]]))
    df_list[[x]] <- rbind(df_list[[x]], union_df[rest_ids, ])
  }
  #rest_ids[1:10]
  filtered_term <- c()

  avgIdx <- list()
  for(TermID in  union_TermID){
    is_sigs <- sapply(df_list, function(x)x[x$ID==TermID,]$p.adjust < 0.05)
    is_true_df <- is_contigous_true_df(is_sigs)
    if(is_true_df$keep){
      filtered_term <- c(filtered_term, TermID)
      avgIdx[[ union_df[union_df$ID==TermID,]$Description ]] <- is_true_df$avgIdx
    }

  }
  if(length(filtered_term) > 5){
    df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_term) )

    nms <- names(df_list)
    df_list <- lapply(names(df_list), function(x) df_list[[x]] %>% mutate(name=x))
    names(df_list) <- nms

    df_list_select <- lapply(1:length(df_list), function(x) df_list[[x]] %>%
                             filter(p.adjust < 0.05) %>%
                             top_n(wt=-log10(p.adjust), n=5) %>%
                             arrange(+log10(p.adjust)))
    df_list_select <- lapply(df_list_select, function(x)x[1:min(5, nrow(x)), ])

    all_names <- as.vector(unlist(sapply(1:length(df_list_select), function(x) (df_list_select[[x]]$ID))))
    pdf_list <- lapply(1:length(df_list), function(x) subset(df_list[[x]], ID %in%all_names))
    mdf <- do.call(rbind, pdf_list)
    pmdf <- mdf[, c("Description", "name", "p.adjust")]
    pmdf$name <- factor(pmdf$name, levels=names(df_list))

    pmtx <- reshape2::dcast(pmdf,  Description ~ name, value.var = "p.adjust")

    rownames(pmtx) <- pmtx$Description
    pmtx$Description <- NULL
    help_mtx <- pmtx
    help_mtx[help_mtx >= 0.05] = 1000
    help_mtx[help_mtx < 0.05] = 1
    help_mtx <- help_mtx[do.call(order, help_mtx),]
    #matrix_list[[pw]] <- pmtx[rownames(help_mtx), ]
    pmtx <- -log10(pmtx)
    pmtx[pmtx>2] = 2
    pmtx <- pmtx[rownames(help_mtx), ]
    ret_mtx <- as.matrix(pmtx)[order(unlist(avgIdx[rownames(help_mtx)])), ]
    return(ret_mtx)
  }
  return(NULL)
}


comb_list <- function(nm){
  m <- combn(nm, 2)
  n <- length(m)/2
  lst <- vector("list", n)

  for (i in 1:n){
    lst[[i]] <- m[1:2, i]
  }
  return(lst)
}


colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}


ggsci_pal <- function(option, ...){
 func_name = glue("pal_{option}")
 func_call = glue('{func_name}(...)')
 assertthat::assert_that(func_name %in% ls("package:ggsci"))
 return(eval(parse(text=func_call)))
}

ggsci_scale_color <- function(option, ...){
 func_name = glue("scale_color_{option}")
 func_call = glue('{func_name}(...)')
 assertthat::assert_that(func_name %in% ls("package:ggsci"))
 return(eval(parse(text=func_call)))
}
ggsci_scale_fill <- function(option, ...){
 func_name = glue("scale_fill_{option}")
 func_call = glue('{func_name}(...)')
 assertthat::assert_that(func_name %in% ls("package:ggsci"))
 return(eval(parse(text=func_call)))
}

URLencode_escape <-function(string){
  string <- URLencode(string, reserved=T)
  string <- stringr::str_replace_all(string, "%", "::::")
  return(string)
}

URLdecode_escape <-function(string){
  string <- stringr::str_replace_all(string, "::::", "%")
  string <- URLdecode(string)
  return(string)
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
  parser <- add_option(parser, c("-j", "--planOfreport"), type="character", default="[\"DEGO_1v1\"]",
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
source(CONFIGFILE)
scrna <- NULL


if (MAKE_ELEMENT == "TRUE"){

  if(identical(cluster,"singleton")){
    cat(paste(date(), blue(" loading: "), red("scrna_phase_singleton.Rds"), "\n"))
    scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
    cat(paste(date(), blue(" loaded: "), red("scrna_phase_singleton.Rds"), "\n"))
  }else{
    cat(paste(date(), blue(" loading: "), red("scrna_phase_comparing.Rds"), "\n"))
    scrna <- readRDS(file=file.path(savedir, "scrna_phase_comparing.Rds"))
    cat(paste(date(), blue(" loaded: "), red("scrna_phase_comparing.Rds"), "\n"))
  }

  scrna$name <- factor(scrna$name, levels=names(data_src))
  scrna$stage <- factor(scrna$stage, levels=unique(stage_lst))
  cluster_viridis_opt = ifelse(
    any(grepl("cluster_color_option",names(viz_conf),fixed = TRUE)),
    viz_conf[["cluster_color_option"]], # Config option
    "d3" # Default
  )

  # replicate_colors
  replicates_viridis_opt = ifelse(
    any(grepl("replicate_color_option",names(viz_conf),fixed = TRUE)),
    viz_conf[["replicate_color_option"]],
    "futurama"
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


  source(glue("{viz_path}/1_quality_report_elements.R"))
  source(glue("{viz_path}/ambientRNA_viz_elements.R"))
  source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  source(glue("{viz_path}/2_batch_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  source(glue("{viz_path}/3_external_markers_elements.R"))
  source(glue("{viz_path}/3_DE_GO-analysis_elements.R"))
  source(glue("{viz_path}/4_DE_GO_1v1_elements.R"))
  source(glue("{viz_path}/4_DE_GO_stageVS_elements.R"))
  source(glue("{viz_path}/4_pathway_1v1_elements.R"))
  source(glue("{viz_path}/4_pathway_stageVS_elements.R"))
  source(glue("{viz_path}/4_Genesets_1v1_elements.R"))
  source(glue("{viz_path}/4_Genesets_stageVS_elements.R"))
  source(glue("{viz_path}/4_progeny_stageVS_elements.R"))

# run necessary generators
  if("QC" %in% EXEC_PLAN) {
    cat(paste(date(), green(" Element: "), red("QC"), "\n"))
    quality_report_elements() ## load scrna innner the function
  }

  if(any(grepl("AmbientRNA",funcs,fixed=TRUE))){
    cat(paste(date(), green(" Element: "), red("AmbientRNA"), "\n"))
    ambientRNA_elements(scrna)
  }
  if("DEs"%in% EXEC_PLAN) {
    cat(paste(date(), green(" Element: "), red("DEs"), "\n"))
    clusters_DEs_elements(scrna)
  }
  if(any(grepl("Clusters_", EXEC_PLAN, fixed=T))){
    cat(paste(date(), green(" Element: "), red("Clusters_integration..."), "\n"))
    batch_clustering_elements(scrna)
  }
  if("Clusters" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("Clusters"), "\n"))
    clustering_elements(scrna)
  }
  if("Singleton" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("Singleton"), "\n"))
    clustering_elements(scrna)
    clusters_DEs_elements(scrna)
  }
  if("EXT_MARKERS" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("EXT_MARKERS"), "\n"))
    external_markers_elements(scrna)
  }
  # FIXME possible problem where all term enrichment analysis is on the same place
  if(length(intersect(c("DEGO","Genesets","progeny","hallmark","KEGG","Reactome"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("DEGO"), "\n"))
    DE_GO_analysis_elements(scrna)
  }
  if("DEGO_1v1" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("DEGO_1v1"), "\n"))
    DEGO_1v1_elements(scrna)
  }
  if("DEGO_stage" %in% EXEC_PLAN){
    cat(paste(date(), green(" Element: "), red("DEGO_stage"), "\n"))
    DEGO_stageVS_elements(scrna)
  }
  if(length(intersect(c("hallmark_1v1","reactome_1v1","kegg_1v1"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("pathway_1v1"), "\n"))
    pathway_1v1_elements(scrna)
  }
  if(length(intersect(c("hallmark_stage","reactome_stage","kegg_stage"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("pathway_stage"), "\n"))
    pathway_stage_elements(scrna)
  }
  if(length(intersect(c("Genesets_stage"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("Genesets_stage"), "\n"))
    Genesets_stageVS_elements(scrna)
  }
  if(length(intersect(c("Genesets_1v1"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("Genesets_1v1"), "\n"))
    Genesets_1v1_elements(scrna)
  }
  if(length(intersect(c("progeny_stage"), EXEC_PLAN) > 0)){
    cat(paste(date(), green(" Element: "), red("progeny_stage"), "\n"))
    progeny_stageVS_elements(scrna)
  }
}

cluster_info <-  build_cluster_info(scrna)

##5. Produce Report
render_func = function(rmd_input_filename, output_filename){
  rmarkdown::render(
    rmd_input_filename,
    output_file=file.path(REPORTDIR,"data",output_filename),
    output_format=c("html_document"),
    quiet=TRUE,
    clean=TRUE,
    params=list(
      scrna=scrna,
      cluster=DEFAULTCLUSTERS,
      cluster_info=cluster_info,
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

dic_Rmd_n_Output <- list(
        "QC"                  =     c(glue("{viz_path}/1_quality_report.Rmd"),        "data_quality"),
        "AmbientRNA"          =     c(glue("{viz_path}/ambientRNA_viz.Rmd"),          "ambient_rna"),
        "DEs"                 =     c(glue("{viz_path}/2_clusters_DEs.Rmd"),          "clusters_DEs"),
        "Clusters"            =     c(glue("{viz_path}/2_clustering.Rmd"),            "clusters"),
        "Singleton"           =     c(glue("{viz_path}/2_clustering.Rmd"),            "clusters"),
        "Clusters_harmony"    =     c(glue("{viz_path}/2_clustering_harmony.Rmd"),    "clusters_harmony"),
        "Clusters_seurat"     =     c(glue("{viz_path}/2_clustering_seurat.Rmd"),     "clusters_seurat"),
        "EXT_MARKERS"         =     c(glue("{viz_path}/3_external_markers.Rmd"),      "external_markers"),
        "DEGO"                =     c(glue("{viz_path}/3_DE_GO-analysis.Rmd"),        "dego"),
        "KEGG"                =     c(glue("{viz_path}/3_KEGG.Rmd"),                  "KEGG"),
        "progeny"             =     c(glue("{viz_path}/3_progeny.Rmd"),               "progeny"),
        "Genesets"            =     c(glue("{viz_path}/3_Genesets.Rmd"),              "Genesets"),
        "hallmark"            =     c(glue("{viz_path}/3_hallmark.Rmd"),              "hallmark"),
        "Reactome"            =     c(glue("{viz_path}/3_Reactome.Rmd"),              "Reactome"),
        "DEGO_stage"          =     c(glue("{viz_path}/4_DE_GO_%s.vs.%s_stageVS.Rmd"),"gv"),
        "DEGO_1v1"            =     c(glue("{viz_path}/4_DE_GO_%s.vs.%s_1v1.Rmd"),    "1vs1"),
        "hallmark_1v1"        =     c(glue("{viz_path}/4_hallmark_1v1.Rmd"),          "hallmark_1vs1"),
        "reactome_1v1"        =     c(glue("{viz_path}/4_reactome_1v1.Rmd"),          "reactome_1vs1"),
        "kegg_1v1"            =     c(glue("{viz_path}/4_kegg_1v1.Rmd"),              "kegg_1vs1"),
        "hallmark_stage"      =     c(glue("{viz_path}/4_hallmark_stageVS.Rmd"),      "hallmark_stageVS"),
        "reactome_stage"      =     c(glue("{viz_path}/4_reactome_stageVS.Rmd"),      "reactome_stageVS"),
        "kegg_stage"          =     c(glue("{viz_path}/4_kegg_stageVS.Rmd"),          "kegg_stageVS"),
        "Genesets_1v1"        =     c(glue("{viz_path}/4_Genesets_1v1.Rmd"),          "Genesets_1vs1"),
        "Genesets_stage"      =     c(glue("{viz_path}/4_Genesets_stageVS.Rmd"),      "Genesets_stageVS"),
        "progeny_stage"      =     c(glue("{viz_path}/4_progeny_stageVS.Rmd"),        "progeny_stageVS"),
        "intUMAPs"            =     c(glue("{viz_path}/interactive_UMAPs.Rmd"),       "interactive_UMAPs")
)




for(i in EXEC_PLAN){
  cat(paste(date(), blue(" Generating: "), red(i), "\n"))
  if(i %ni% names(dic_Rmd_n_Output)){
    message("viz ", i, " is not implemented!!")
    next
  }
  rmd_n_output <- dic_Rmd_n_Output[[i]]
  rmd <- rmd_n_output[1]
  output <- rmd_n_output[2]
  if(output == "1vs1"){
    for(apair in comb_list(names(data_src))){
      render_func(sprintf(rmd, apair[1], apair[2]), glue("{output}_{apair[1]}.vs.{apair[2]}.html"))
    }
  }else if(output == "gv"){
    for(apair in comb_list(unique(stage_lst))){
      render_func(sprintf(rmd, apair[1], apair[2]), glue("{output}_{apair[1]}.vs.{apair[2]}.html"))
    }
  }else{
    render_func(rmd, output)
  }

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

