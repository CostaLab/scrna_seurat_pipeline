is_contigous_true_df <- function(is_sigs){
  ret_df <- data.frame(keep=FALSE, avgIdx=-1)
  if(any(is_sigs) & table(is_sigs)['TRUE'] == 1){
        ret_df$keep=TRUE
        ret_df$avgIdx = which(is_sigs == TRUE)
        return(ret_df)
  }
  return(ret_df)
}

scrna <- readRDS(file.path(savedir, "scrna_phase_comparing.Rds"))
de_cluster_name <- paste0("de_", cluster)
go_cluster_name <- paste0("go_", cluster)

if(de_cluster_name %ni% names(scrna@tools)){
  stop(glue("ERROR:DE hasn't been calculated for cluster:{cluster}\n Please run [scrna_markergenes]!!!"))
}
if(go_cluster_name %ni% names(scrna@tools)){
  stop(glue("ERROR:GO hasn't been calculated for cluster:{cluster}\n Please run [scrna_go]!!!"))
}
if(any(grepl("hallmark",funcs,fixed=TRUE))){
  hallmark_cluster_name <- paste0("hallmark_", cluster)
  if(hallmark_cluster_name %ni% names(scrna@tools)){
    stop(glue("ERROR:hallmark hasn't been calculated for cluster {cluster}\n Please run [scrna_hallmark]!!!"))
  }
}
if(any(grepl("KEGG",funcs,fixed=TRUE))){
  kegg_cluster_name <- paste0("kegg_", cluster)
  if(kegg_cluster_name %ni% names(scrna@tools)){
    stop(glue("ERROR:kegg hasn't been calculated for cluster {cluster}\n Please run [scrna_kegg]!!!"))
  }
}
if(any(grepl("Reactome",funcs,fixed=TRUE))){
  reactome_cluster_name <- paste0("reactome_", cluster)
  if(reactome_cluster_name %ni% names(scrna@tools)){
    stop(glue("ERROR:reactome hasn't been calculated for cluster {cluster}\n Please run [scrna_reactome]!!!"))
  }
}

# DE plots
cluster_de <- scrna@tools[[de_cluster_name]]
cluster_de <- cluster_de[sapply(cluster_de, function(m) nrow(m) >0)]

cluster_de_top10 <- lapply(cluster_de, function(x) {
    x %>% top_n(10, avg_logFC) %>% arrange(-avg_logFC)
})


## top10 DE heatmaps
genes <- as.vector(unlist(sapply(cluster_de_top10, function(x)x$gene)))
scrna <- ScaleData(scrna, rownames(scrna))
help_sort_func <- ifelse(
  all.is.numeric(unique(scrna@meta.data[, cluster])),
  function(x) as.numeric(x)-1,
  as.character
)

scrna@meta.data[, cluster] <- help_sort_func(scrna@meta.data[, cluster])

plthm = DoHeatmap(
  scrna,
  features=genes,
  group.by = cluster,
  disp.min = -2,
  disp.max = 2,
  slot = "scale.data",
  assay = "RNA",
  raster = FALSE,
  combine= T
) +
ggtitle("Marker genes for each cluster") +
NoLegend()

save_ggplot_formats(
  plt=plthm,
  base_plot_dir=report_plots_folder,
  plt_name=paste0("heatmap_top10-de-genes_cluster-",cluster),
  width=13, height=12
)


# GeneBarPlots
## Plot the top 10 DE genes in each cluster.
plots = list()

help_sort_func <- ifelse(all.is.numeric(names(cluster_de)), as.numeric, function(x){x})

for (id in sort(help_sort_func(names(cluster_de)))) {
  id = as.character(id)
  cluster_genes <- cluster_de_top10[[id]]
  x_lim = max(abs(cluster_de[[id]]$avg_logFC))
  x_lim <- c(-x_lim, x_lim)
  plots[[id]] <- GeneBarPlot(cluster_de[[id]], xlim = x_lim, main = id)
}

if(length(plots) > 0){
	for (i in seq(1, length(plots), by=4)){
	  ni = min(i+3, length(plots))
	  plt <-plot_grid(plotlist=plots[i:ni], ncol=4)
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("top10-de-genes_per_cluster-",cluster,"_p",i,"-",ni),
      width=13, height=12
    )
	}
}


# volcano plots
help_sort_func <- ifelse(all.is.numeric(names(cluster_de)), as.numeric, function(x){x})

for (id in sort(help_sort_func(names(cluster_de)))) {
  id = as.character(id)
  a_de <- cluster_de[[id]]
  a_de$log2FC <- a_de$avg_logFC / log(2)
  up <- nrow(a_de %>% filter(log2FC>= 1 & p_val_adj<=0.05) )
  down <- nrow(a_de %>% filter(log2FC <= -1 & p_val_adj<=0.05))
  plt <- EnhancedVolcano(
    a_de,
    x="log2FC",
    y = "p_val_adj",
    lab=rownames(a_de),
    pointSize = 1.0,
    pCutoff = 0.05,
    title=glue("Volcano {id}"),
    subtitle=glue("up:{up} down:{down}")
  )
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("volcanoplot_deg_cluster-",cluster,"_id-",id),
    width=8, height=8
  )
}

## DE genes on UMAP plot
Idents(scrna) <- cluster
for (i in names(cluster_de) ){
  #plots = list()
  print(sprintf("Cluster %s:", i))
  ps <- FeaturePlot(
    scrna,
    features = cluster_de_top10[[as.character(i)]]$gene,
    label=T,
    label.size=2,
    cols = c("lightgrey", "red"),
    reduction = "DEFAULT_UMAP"
  )
  save_ggplot_formats(
    plt=ps,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_featureplot_top10deg_cluster-",cluster,"_id-",i),
    width=13, height=10
  )
}

#############################################################
## Term enrichment analysis (GO, hallmark, KEGG, Reactome) ##
#############################################################
enrch_analysis_vec = go_cluster_name
if(any(grepl("hallmark",funcs,fixed=TRUE))) enrch_analysis_vec = c(enrch_analysis_vec, hallmark_cluster_name)
if(any(grepl("KEGG",funcs,fixed=TRUE))) enrch_analysis_vec = c(enrch_analysis_vec, kegg_cluster_name)
if(any(grepl("Reactome",funcs,fixed=TRUE))) enrch_analysis_vec = c(enrch_analysis_vec, reactome_cluster_name)


for(enrich_cluster_name in enrch_analysis_vec){



  enrich_type = gsub(paste0("_",cluster,"$"),"",enrich_cluster_name)
  ## TERM up analysis
  term_up_list <- if(enrich_cluster_name == go_cluster_name){scrna@tools[[go_cluster_name]]$goup}
    else if(enrich_cluster_name == hallmark_cluster_name){scrna@tools[[hallmark_cluster_name]]$hallmarkup}
    else if(enrich_cluster_name == kegg_cluster_name){scrna@tools[[kegg_cluster_name]]$keggup}
    else if(enrich_cluster_name == reactome_cluster_name){scrna@tools[[reactome_cluster_name]]$reactomeup}


  ## TERM down analysis
  term_down_list <- if(enrich_cluster_name == go_cluster_name){scrna@tools[[go_cluster_name]]$godown}
    else if(enrich_cluster_name == hallmark_cluster_name){scrna@tools[[hallmark_cluster_name]]$hallmarkdown}
    else if(enrich_cluster_name == kegg_cluster_name){scrna@tools[[kegg_cluster_name]]$keggdown}
    else if(enrich_cluster_name == reactome_cluster_name){scrna@tools[[reactome_cluster_name]]$reactomedown}

  for(term_direction in c("up","down")){

    term_direction_list = switch(term_direction, "up"=term_up_list, "down"=term_down_list)

    df_list <- lapply(1:length(term_direction_list), function(x) term_direction_list[[x]]@result)
    names(df_list) <- names(term_direction_list)
    for(i in 1:length(df_list)){
        assertthat::assert_that( all(rownames(df_list[[i]]) == df_list[[i]]$ID))
    }
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

      pmtx <- reshape2::dcast(pmdf,  Description ~ name)

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
      col_fun <-  circlize::colorRamp2(c(0, 0.5, 2), c("purple", "black", "yellow"))
      plthm <- Heatmap(as.matrix(pmtx)[order(unlist(avgIdx[rownames(help_mtx)])), ],
             name = glue("-log10(padjust)"),
             cluster_columns = F,
             cluster_rows = F,
             show_row_names=T,
             col=col_fun,
             heatmap_legend_param = list(direction = "horizontal")
      )

      png(
        filename=file.path(report_plots_folder_png,paste0("heatmap_",enrich_type,"-",term_direction,"_cluster-",cluster,".png")),
        width=13,
        height=20,
        type="cairo",
        units="in",
        res=300
      )
      draw(plthm, heatmap_legend_side = "top")
      dev.off()
      pdf(
        file=file.path(report_plots_folder_pdf,paste0("heatmap_",enrich_type,"-",term_direction,"_cluster-",cluster,".pdf")),
        width=13,
        height=20
      )
      draw(plthm, heatmap_legend_side = "top")
      dev.off()
    }

    ### direction genes top 10
    term_plot_list = lapply(
      term_direction_list,
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
        ggtitle(paste(enrich_type,term_direction,",", names(y)[i])) +
        theme_minimal() +
        theme(axis.text.y  = element_text(size=20), axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
        labs(y = "") +
        scale_y_continuous(name="-log10(p-value)") +
        scale_x_discrete(breaks = col[,2], labels = col[,3]) +
        #scale_x_discrete(name= "") +
        coord_flip()
      },
      y=term_plot_list
    )

    #plt = plot_grid(plotlist=plots, ncol=2)
    if(length(plots) > 0){
        for (i in seq(1, length(plots), by=2)){
          ni = min(i+1, length(plots))
          plt <-plot_grid(plotlist=plots[i:ni], ncol=2)

          save_ggplot_formats(
            plt=plt,
            base_plot_dir=report_plots_folder,
            plt_name= paste0(enrich_type,"_",term_direction,"_genes_barplot_cluster-",cluster, "_p",i, "-", ni),
            width=16, height=10)
        }
    }
  }
}
