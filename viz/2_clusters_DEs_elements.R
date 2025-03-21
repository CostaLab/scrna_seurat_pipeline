clusters_DEs_elements <- function(scrna){
  col_string=""
  if(identical(cluster,"singleton")){
    col_string = "RNA_snn_res."
  }else{
    col_string = "integrated_snn_res."
    if(INTEGRATION_OPTION=="harmony"){
      col_string = "RNA_snn_res."
    }
  }
  cluster_de_list  <- seutools_partition(scrna,
                                         partition="de_batch",
                                         save_dir=SAVE_DIR,
                                         allinone=ALLINONE)

  ## Plot the top 10 DE genes in each cluster.
  # FIXME shouldnt it be vector that the user defined?
  names(cluster_de_list) <- as.character(seq(0.1, 0.8, 0.1)) ##------
  for (resolution in seq(0.1, 0.8, 0.1)){

    cluster_de <- cluster_de_list[[as.character(resolution)]]
    cluster_de <- cluster_de[sapply(cluster_de, function(m) nrow(m) >0)]
    cluster_de_top10 <- lapply(cluster_de, function(x) {
        if("avg_logFC" %in% names(x)){ ## compatible for seurat3
          x$avg_log2FC <- x$avg_logFC/log(2)
        }
        x %>% top_n(10, avg_log2FC) %>% arrange(-avg_log2FC)
    })

    plots = list()

    help_sort_func <- ifelse(
      all.is.numeric(names(cluster_de)),
      as.numeric,
      function(x){x}
    )

    for (id in sort(help_sort_func(names(cluster_de)))) {
      id = as.character(id)
      cluster_genes = cluster_de_top10[[id]]
      if("avg_logFC" %in% names(cluster_de[[id]])){## compatible for seurat3
        cluster_de[[id]]$avg_log2FC <-cluster_de[[id]]$avg_logFC / log(2)
      }
      x_lim <- max(abs(cluster_de[[id]]$avg_log2FC))
      x_lim <- c(-x_lim, x_lim)
      plots[[id]] <- GeneBarPlot(
        cluster_de[[id]],
        xlim = x_lim,
        main = paste0("Cluster: ", id)
      )
    }

    ps <- plot_grid(plotlist = plots)
    plt_title <- ggdraw() + 
      draw_label(sprintf("resolution %s", resolution), fontface = "bold")

    plt <- plot_grid(plt_title, ps, ncol = 1, rel_heights = c(0.1, 1))

    save_ggplot_formats(
      plt = plt,
      base_plot_dir = report_plots_folder,
      plt_name = paste0("clusterDE_genes_resolution-", resolution),
      width = 13, height = 20
    )
  }


  ## DE genes on heatmap
  DefaultAssay(scrna) <- "RNA"
  scrna <- Seurat::ScaleData(scrna,  rownames(scrna))

  cluster_de_list  <- seutools_partition(scrna,
                                         partition="de_batch",
                                         save_dir=SAVE_DIR,
                                         allinone=ALLINONE)


  names(cluster_de_list) <- as.character(seq(0.1, 0.8, 0.1))

  for (resolution in seq(0.1, 0.8, 0.1)){
    #cluster.de <- cluster.de.list[[as.character(resolution)]]
    cluster_de_top8 <- lapply(cluster_de, function(x) {
        if("avg_logFC" %in% names(x)){ ## compatible for seurat3
          x$avg_log2FC <- x$avg_logFC / log(2)
        }
        x %>% top_n(8, avg_log2FC) %>% arrange(-avg_log2FC)
    })

    cluster_de_top8_combine <- do.call(rbind, cluster_de_top8)
    genes <- unique(cluster_de_top8_combine$gene)

    col_def <- rev(
      ggsci_pal(option = cluster_viridis_opt)(
        length(unique(scrna@meta.data[,sprintf(paste0(col_string,"%.1f"), resolution)]))
      )
    )

    scrna_to_plot <- function(scrna){
     if(ncol(scrna)<=20000){
       return(scrna)
     }
     return(subset(scrna, cells=sample(colnames(scrna))[1:20000]))
    }


    plt <- DoHeatmap(
      ## >30k will failed to plot, here we subset when cell number > 20,000
      object = scrna_to_plot(scrna),
      features = genes,
      group.by = sprintf(paste0(col_string, "%.1f"), resolution),
      group.colors = col_def,
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
      plt = plt,
      base_plot_dir = report_plots_folder,
      plt_name = paste0("heatmapDE_genes_resolution-", resolution),
      width = 10, height = 13
    )
  }
}

