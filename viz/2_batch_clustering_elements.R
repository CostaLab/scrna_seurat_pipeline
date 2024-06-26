
batch_clustering_elements <- function(scrna){
  ## Clusters Resolution
  #if(identical(cluster,"singleton")){
  #  scrna <- load_object(file_name = file.path(savedir, "scrna_phase_singleton.Rds"))
  #}else{
  #  scrna <- load_object(file_name = file.path(savedir, "scrna_phase_comparing.Rds"))
  #}

  for(cluster_use in available_clusters){

    message(paste0("### Producing elements for cluster: ",cluster_use))

    if(cluster_use %ni% names(scrna@meta.data)){
      stop(glue("ERROR:There's no this {cluster_use} slot, please check!!!"))
    }

    pref_def = "integrated_snn_res."
    if(cluster_use == "harmony_inte_clusters") pref_def = "RNA_snn_res."
    if(cluster_use == "singleton") pref_def = "RNA_snn_res."

    umap_reduction = "DEFAULT_UMAP"
    if(cluster_use == "harmony_inte_clusters") umap_reduction = "harmony_UMAP"
    if(cluster_use == "seurat_inte_clusters") umap_reduction = "INTE_UMAP"
    if(cluster_use == "singleton") umap_reduction = "SINGLE_UMAP"

    message("### Making cluster tree")
    plt = clustree(scrna, prefix = pref_def)
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("clustree_resolution_",cluster_use),
      width=13, height=10
    )

    message("### Making umap resolution list")
    # UMAP resolution list
    # FIXME should be the resolution vector the user defined right?
    nms <- as.character(seq(0.1, 0.8, 0.1))
    plist <- list()
    for(nm in nms){

      col_def <- rev(
        ggsci_pal(option = cluster_viridis_opt)(
          length(unique(scrna@meta.data[,paste0(pref_def, nm)]))
        )
      )

      plist[[nm]] <- DimPlot(
        scrna,
        reduction = umap_reduction,
        group.by =  paste0(pref_def, nm),
        label=TRUE,
        label.size=8,
        cols=col_def
      ) + ggtitle(sprintf("resolution %s", nm))

    }

    plt = patchwork::wrap_plots(plist, ncol=2)
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_resolution_list_",cluster_use),
      width=15, height=20
    )
  }
}

