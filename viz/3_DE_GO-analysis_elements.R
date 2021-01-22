

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
    reduction = "INTE_UMAP"
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

    term_direction_list = ifelse(term_direction == "up", term_up_list, term_down_list)

    df_list <- lapply(names(term_direction_list), function(x) term_direction_list[[x]]@result)
    names(df_list) <- names(term_direction_list)
    intersect_termID <- Reduce(intersect, lapply(df_list, function(x) x$ID))
    filtered_term <- c()
    for(termid in  intersect_termID){
        is_sigs <- sapply(df_list, function(x)x[x$ID==termid,]$p.adjust < 0.05)
        if(any(is_sigs) & table(is_sigs)["TRUE"] == 1){
            filtered_term <- c(filtered_term, termid)
        }
    }

    if(length(filtered_term) > 10){
      df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_term) )
      df_list <- lapply(names(df_list), function(x) df_list[[x]] %>% mutate(name=x))
      mdf <- do.call(rbind, df_list)
      pmdf <- mdf[, c("Description", "name", "p.adjust")]

      pmtx <- reshape2::dcast(pmdf,  Description ~ name)

      rownames(pmtx) <- pmtx$Description
      pmtx$Description <- NULL
      help_mtx <- pmtx
      help_mtx[help_mtx >= 0.05] = 1000
      help_mtx[help_mtx < 0.05] = 1
      help_mtx <- help_mtx[do.call(order, help_mtx),]
      pmtx <- -log10(pmtx)
      pmtx[pmtx>2] = 2
      pmtx <- pmtx[rownames(help_mtx), ]
      col_fun <-  circlize::colorRamp2(c(0, 1, +2), c("purple", "black", "yellow"))
      plthm <- Heatmap(
        as.matrix(pmtx),
        name = "-log10(padjust)",
        cluster_columns = F,
        cluster_rows = F,
        show_row_names=T,
        col=col_fun
      )
      png(
        filename=file.path(report_plots_folder_png,paste0("heatmap_",enrich_type,"-",term_direction,"_cluster-",cluster,".png")),
        width=13,
        height=20,
        type="cairo",
        units="in",
        res=300
      )
      draw(plthm, heatmap_legend_side = "left")
      dev.off()
      pdf(
        file=file.path(report_plots_folder_pdf,paste0("heatmap_",enrich_type,"-",term_direction,"_cluster-",cluster,".pdf")),
        width=13,
        height=20
      )
      draw(plthm, heatmap_legend_side = "left")
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
          if(str_length(y) > 50){
            hs <- digest(y, "crc32")
            y = paste(substr(y, 1, 40), hs)}
            return(y) }))
        #log_m$names = df$Description
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
        geom_bar(stat="identity", fill="#3399CC", color="grey50") +
        ggtitle(paste(enrich_type,term_direction,",", names(y)[i])) +
        theme(axis.text.y  = element_text(size=20)) +
        scale_y_continuous(name="-log10(p-value)") +
        scale_x_discrete(name= "") +
        coord_flip()
      },
      y=term_plot_list
    )

    plt = plot_grid(plotlist=plots, ncol=2)
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0(enrich_type,"_",term_direction,"_genes_barplot_cluster-",cluster),
      width=12, height=20
    )
  }
}
