
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

run_shell <- function(cmd){
  system(cmd)
}

save_ggplot_formats = function(
  plt, base_plot_dir, plt_name, create_plot_subdir=TRUE,
  formats=c("png", "pdf"), units="in", width=20, height=20,
  type="cairo", res = 300,
  plot_obj="ggplot",
  ...
){
  # if theres a plot and basedir
  if(!is.null(base_plot_dir) & !is.null(plt)){
    # for each format
    for(fmt in formats){
      f_path_fmt = file.path(base_plot_dir, paste0(plt_name, ".", fmt))
      if(create_plot_subdir) dir.create(file.path(base_plot_dir,fmt),recursive=TRUE,showWarnings=FALSE)
      if(dir.exists(file.path(base_plot_dir,fmt))) f_path_fmt = file.path(base_plot_dir,fmt,paste0(plt_name,".",fmt))
      if(plot_obj == "ggplot"){
        if(fmt == "png"){
          ggplot2::ggsave(
            filename = f_path_fmt,
            plot = plt,
            device = fmt,
            units = units,
            width = width,
            height = height,
            type = type, ...
          )
        }else{
          ggplot2::ggsave(
            filename = f_path_fmt,
            plot = plt,
            device = fmt,
            units = units,
            width = width,
            height = height, ...
          )
        }
      }else
        if(fmt == "png"){
          png(
            filename = f_path_fmt,
            units = units,
            width = width,
            height = height,
            type = type,
            res = res
          )
          draw(plt, ...)
          dev.off()
        }else{
          pdf(
            file = f_path_fmt,
            width = width,
            height = height
          )
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
    de.data$avg_log2FC <- de.data$avg_logFC / log(2)
  }
  if (any(colnames(de.data) == "cluster")) {
    top5.up <-
      de.data %>%
      group_by(cluster) %>%
      top_n(10, avg_log2FC) %>%
      filter(avg_log2FC > 0) %>%
      arrange(-avg_log2FC)
    top5.dn <-
      de.data %>%
      group_by(cluster) %>%
      top_n(10, -avg_log2FC) %>%
      filter(avg_log2FC < 0) %>%
      arrange(-avg_log2FC)
  } else {
    top5.up <-
      de.data %>%
      top_n(10, avg_log2FC) %>%
      filter(avg_log2FC > 0) %>%
      arrange(-avg_log2FC)
    top5.dn <-
      de.data %>%
      top_n(10, -avg_log2FC) %>%
      filter(avg_log2FC < 0) %>%
      arrange(-avg_log2FC)
  }
  
  top.up.dn <- rbind(top5.up, top5.dn)
  top.up.dn$gene <- make.unique(top.up.dn$gene)
  top.up.dn$type <- ifelse(top.up.dn$avg_log2FC > 0, "positive", "negative")
  top.up.dn$type <- factor(top.up.dn$type, levels = c("positive", "negative"))
  
  g <-
    ggplot(
      data = top.up.dn,
      aes(x = gene, y = avg_log2FC, fill = type)
    ) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = rev(top.up.dn$gene)) +
    theme_minimal() +
    theme(legend.position = "none", axis.text = element_text(size = 15)) +
    scale_fill_manual(values = c(positive = pos_color, negative = neg_color)) +
    coord_flip()
  if (!is.null(main)) {
    g <- g + ggtitle(paste0(main))
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
    p <- ggplot(data = counts, aes(x = x, y = Cells, fill = fill)) +
      geom_bar(stat = "identity",position = position_dodge())
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
      if (is.list(df) && length(df) == 0){
        log_m = as.data.frame(list())
        return(log_m)
      }
      log_m = as.data.frame(-log10(df$p.adjust))
      log_m$names = as.factor(
        sapply(
          df$Description,
          function(y){
            y <- as.character(trimws(y))
            return(y)
          }
        )
      )
      log_m$show_names = as.factor(
        sapply(
          df$Description,
          function(y){
            y <- as.character(trimws(y))
            y <- ifelse(nchar(y) <= 33,  y, paste0(substr(y, 1, 30), "..."))
            return(y)
          }
        )
      )
      log_m <- log_m[order(log_m[,1],decreasing = TRUE),]
      showCatetermry <- min(length(log_m[,1]), 10)
      log_m <- log_m[1:showCatetermry, ]
      log_m <- log_m[order(log_m[,1],decreasing = FALSE),]
      return(log_m)
    }
  )

  ### direction genes plot
  plots <- lapply(
    seq_along(term_plot_list),
    function(y, i) {
      term_df <- y[[i]]
      if(length(term_df) == 0) return(NULL)

      ggplot(
        term_df,
        aes(x = reorder(x = term_df[,2], term_df[,1]), y = term_df[,1])
      ) +
      geom_bar(
        stat = "identity", fill = pos_color
      ) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = 2,
        size = 1,
        color = "grey50"
      ) +
      ggtitle(paste0("Cluster: ", names(y)[i])) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.y = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      labs(y = "") +
      scale_y_continuous(name = "-log10(p-value)") +
      scale_x_discrete(breaks = term_df[,2], labels = term_df[,3]) +
      coord_flip()
    },
    y = term_plot_list
  )
  plots <- Filter(Negate(is.null), plots)
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

    df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_term))

    nms <- names(df_list)
    df_list <- lapply(names(df_list), function(x) df_list[[x]] %>% mutate(name=x))
    names(df_list) <- nms

    df_list_select <- lapply(
      1:length(df_list),
      function(x){
        df_list[[x]] %>%
          filter(p.adjust < 0.05) %>%
          top_n(wt = -log10(p.adjust), n = 5) %>%
          arrange(log10(p.adjust))
      }
    )
    df_list_select <- lapply(df_list_select, function(x)x[1:min(5, nrow(x)), ])

    all_names <- as.vector(unlist(sapply(1:length(df_list_select), function(x) (df_list_select[[x]]$ID))))
    pdf_list <- lapply(1:length(df_list), function(x) subset(df_list[[x]], ID %in%all_names))
    mdf <- do.call(rbind, pdf_list)
    pmdf <- mdf[, c("Description", "name", "p.adjust")]
    pmdf$name <- factor(pmdf$name, levels = names(df_list))

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


