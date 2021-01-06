

scrna <- readRDS(file.path(savedir, "scrna_phase_comparing.Rds"))
de_cluster_name <- paste0("de_", cluster)
go_cluster_name <- paste0("go_", cluster)

if(de_cluster_name %ni% names(scrna@tools) | go_cluster_name %ni% names(scrna@tools)){
  stop(glue("ERROR:DE&GO hasn't been calculated for cluster:{cluster}\n Please run [scrna_markergenes, scrna_go]!!!"))
}

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






## GO up analysis
go_up_list <-scrna@tools[[go_cluster_name]]$goup

df_list <- lapply(names(go_up_list), function(x) go_up_list[[x]]@result)
names(df_list) <- names(go_up_list)
intersect_GoID <- Reduce(intersect, lapply(df_list, function(x) x$ID))
filtered_go <- c()
for(goid in  intersect_GoID){
    is_sigs <- sapply(df_list, function(x)x[x$ID==goid,]$p.adjust < 0.05)
    if(any(is_sigs) & table(is_sigs)["TRUE"] == 1){
        filtered_go <- c(filtered_go, goid)
    }
}

if(length(filtered_go) > 10){
  df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_go) )
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
    filename=file.path(report_plots_folder_png,paste0("heatmap_go-up_cluster-",cluster,".png")),
    width=13,
    height=20,
    units="in",
    res=300
  )
  draw(plthm, heatmap_legend_side = "left")
  dev.off()
  pdf(
    filename=file.path(report_plots_folder_pdf,paste0("heatmap_go-up_cluster-",cluster,".pdf")),
    width=13,
    height=20
  )
  draw(plthm, heatmap_legend_side = "left")
  dev.off()
}



### up genes top 10
go_plot_list = lapply(go_up_list, function(x){
  df = x@result
  if (is.list(df) && length(df)==0){
    log_m = as.data.frame(list())
    return(log_m)
  }
  log_m = as.data.frame(-log10(df$p.adjust))
  log_m$names =as.factor(sapply(df$Description, function(y){
    y <- as.character(trimws(y))
    if(str_length(y) > 50){
      hs <- digest(y, "crc32")
      y = paste(substr(y, 1, 40), hs)}
      return(y) }))
  #log_m$names = df$Description
  log_m <- log_m[order(log_m[,1],decreasing = TRUE),]
  showCategory = min(length(log_m[,1]), 10)
  log_m <- log_m[1:showCategory, ]
  log_m <- log_m[order(log_m[,1],decreasing = FALSE),]
  return(log_m)
  }
)

### up genes plot
plots <- lapply(
  seq_along(go_plot_list),
  function(y, i) {
    col <- y[[i]]
    if(length(col) == 0)
      return(NULL)
    ggplot(col, aes(reorder(x=col[,2], col[,1]), y=col[,1])) +
    geom_bar(stat="identity", fill="#3399CC", color="grey50") +
    ggtitle(paste("GO Up, ", names(y)[i])) +
    theme(axis.text.y  = element_text(size=20)) +
    scale_y_continuous(name="-log10(p-value)") +
    scale_x_discrete(name= "") +
    coord_flip()
  },
  y=go_plot_list
)

plt = plot_grid(plotlist=plots, ncol=2)
save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name=paste0("go_up_genes_barplot_cluster-",cluster),
  width=12, height=20
)





## GO down analysis
go_down_list <-scrna@tools[[go_cluster_name]]$godown

df_list <- lapply(names(go_down_list), function(x) go_down_list[[x]]@result)
names(df_list) <- names(go_down_list)
intersect_GoID <- Reduce(intersect, lapply(df_list, function(x) x$ID))
filtered_go <- c()
for(goid in  intersect_GoID){
    is_sigs <- sapply(df_list, function(x)x[x$ID==goid,]$p.adjust < 0.05)
    if(any(is_sigs) & table(is_sigs)["TRUE"] == 1){
        filtered_go <- c(filtered_go, goid)
    }
}

if(length(filtered_go) > 10){
  df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_go) )
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
    filename=file.path(report_plots_folder_png,paste0("heatmap_go-down_cluster-",cluster,".png")),
    width=13,
    height=20,
    units="in",
    res=300
  )
  draw(plthm, heatmap_legend_side = "left")
  dev.off()
  pdf(
    filename=file.path(report_plots_folder_pdf,paste0("heatmap_go-down_cluster-",cluster,".pdf")),
    width=13,
    height=20
  )
  draw(plthm, heatmap_legend_side = "left")
  dev.off()
}


### down genes top 10
go_plot_list = lapply(
  go_down_list,
  function(x){
    #df = fortify(x, showCategory=Inf)
    df = x@result
    if (is.list(df) && length(df)==0){
      log_m = as.data.frame(list())
      return(log_m)
    }
    log_m = as.data.frame(-log10(df$p.adjust))
    log_m$names =as.factor(sapply(df$Description, function(y){
      y <- as.character(trimws(y))
      if(str_length(y) > 50){
        hs <- digest(y, "crc32")
        y = paste(substr(y, 1, 40), hs)}
        return(y) }))
    #log_m$names = df$Description
    log_m <- log_m[order(log_m[,1],decreasing = TRUE),]
    showCategory = min(length(log_m[,1]), 10)
    log_m <- log_m[1:showCategory, ]
    log_m <- log_m[order(log_m[,1],decreasing = FALSE),]
    return(log_m)
  }
)

### down genes plot
plots<-lapply(
  seq_along(go_plot_list),
  function(y, i) {
    col <- y[[i]]
    if(length(col) == 0)
      return(NULL)
    ggplot(col, aes(reorder(x=col[,2], col[,1]), y=col[,1])) +
    geom_bar(stat="identity", fill= "#3399CC", color="grey50") +
    ggtitle(paste("GO Down", names(y)[i])) +
    theme(axis.text.y  = element_text(size=20)) +
    scale_y_continuous(name="-log10(p-value)") +
    scale_x_discrete(name= "") +
    coord_flip()
  },
  y=go_plot_list
)

plt = plot_grid(plotlist=plots, ncol=2)
save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name=paste0("go_down_genes_barplot_cluster-",cluster),
  width=12, height=20
)



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
