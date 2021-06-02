
col_string=""
if(identical(cluster,"singleton")){
  # scrna <- readRDS(file=file.path(savedir, "scrna_singleton_markergenes.Rds"))
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
  col_string = "RNA_snn_res."
}else{
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_comparing.Rds"))
  col_string = "integrated_snn_res."
}
cluster_de_list <- scrna@tools$de_batch

## Plot the top 10 DE genes in each cluster.
# FIXME shouldnt it be vector that the user defined?
names(cluster_de_list) <- as.character(seq(0.1, 0.8, 0.1)) ##------
for (resolution in seq(0.1, 0.8, 0.1)){
  cluster_de <- cluster_de_list[[as.character(resolution)]]
  cluster_de <- cluster_de[sapply(cluster_de, function(m) nrow(m) >0)]
  cluster_de_top10 <- lapply(cluster_de, function(x) {
      x %>% top_n(10, avg_logFC) %>% arrange(-avg_logFC)
  })

  plots = list()
  help_sort_func <- ifelse(all.is.numeric(names(cluster_de)), as.numeric, function(x){x})
  for (id in sort(help_sort_func(names(cluster_de)))) {
    id = as.character(id)
    cluster_genes = cluster_de_top10[[id]]
    x_lim = max(abs(cluster_de[[id]]$avg_logFC))
    x_lim = c(-x_lim, x_lim)
    plots[[id]] = GeneBarPlot(cluster_de[[id]], xlim = x_lim, main = id)
  }

  ps <- plot_grid(plotlist=plots)
  title <- ggdraw() + draw_label(sprintf("resolution %s", resolution), fontface='bold')

  plt <- plot_grid(title, ps, ncol=1, rel_heights=c(0.1, 1))
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("clusterDE_genes_resolution-", resolution),
    width=13, height=20
  )
}


## DE genes on heatmap
DefaultAssay(scrna) <- "RNA"
scrna <- Seurat::ScaleData(scrna,  rownames(scrna))

cluster_de_list <- scrna@tools$de_batch
names(cluster_de_list) <- as.character(seq(0.1, 0.8, 0.1))

for (resolution in seq(0.1, 0.8, 0.1)){
  #cluster.de <- cluster.de.list[[as.character(resolution)]]
  cluster_de_top8 <- lapply(cluster_de, function(x) {
      x %>% top_n(8, avg_logFC) %>% arrange(-avg_logFC)
  })

  cluster_de_top8_combine <- do.call(rbind, cluster_de_top8)
  genes <- unique(cluster_de_top8_combine$gene)

  col_def <- rev(
    viridis::viridis_pal(option = cluster_viridis_opt)(
      length(unique(scrna@meta.data[,sprintf(paste0(col_string,"%.1f"), resolution)]))
    )
  )

  plt <- DoHeatmap(
    object = scrna,
    features = genes,
    group.by = sprintf(paste0(col_string,"%.1f"), resolution),
    group.colors=col_def,
    disp.min = -2,
    disp.max = 2,
    slot = "scale.data",
    assay = "RNA",
    raster = FALSE,
    combine = TRUE
  ) +
  ggtitle(sprintf("resolution: %.1f", resolution)) +
  NoLegend()

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("heatmapDE_genes_resolution-", resolution),
    width=10, height=13
  )
}
