####################################################
# doublets 
####################################################
doubletdetection_viz_elements <- function(scrna){
  # We have three options.
  # If doublet_switch == "display", we have all the cells still in the Seurat object.
  # We plot using DEFAULT_UMAP when available, otherwise fall back to harmony/INTE UMAP.
  # If doublet_switch == "on", we have removed the doublets from the Seurat object.
  # We can use the output from the doublet detection function (scrna_DoubletAnnotated.Rds) and plot the UMAP we created for this purpose.
  # If doublet_switch == "off", we just put out a message stating that doublet detection was not performed.
  if(doublet_switch == "display"){
    scrna <- load_object(file_name = file.path(savedir, "scrna_phase_comparing.Rds"))
    scrna <- safe_join_layers(scrna)  # v5: join layers
    if ("Doublet_classifications" %in% names(scrna@meta.data)) {
      reduction_use <- "DEFAULT_UMAP"
      if (!("DEFAULT_UMAP" %in% names(scrna@reductions))) {
        if ("harmony_UMAP" %in% names(scrna@reductions)) {
          reduction_use <- "harmony_UMAP"
        } else if ("INTE_UMAP" %in% names(scrna@reductions)) {
          reduction_use <- "INTE_UMAP"
        }
      }
      if (reduction_use %in% names(scrna@reductions)) {
        plt <- DimPlot(scrna, group.by = "Doublet_classifications", reduction = reduction_use)
        save_ggplot_formats(
          plt=plt,
          base_plot_dir=report_plots_folder,
          plt_name="doublets_umap",
          width=9, height=7
        )
      } else {
        message("No suitable UMAP reduction found for doublet display plot; skipping.")
      }
    }
  } else if(doublet_switch == "on"){
    scrna <- load_object(file_name = file.path(savedir, "scrna_DoubletAnnotated.Rds"))
    scrna <- safe_join_layers(scrna)  # v5: join layers
    if ("Doublet_classifications" %in% names(scrna@meta.data)) {
      plt <- DimPlot(scrna, group.by = "Doublet_classifications", reduction = "DOUBLET_UMAP")
      save_ggplot_formats(
        plt=plt,
        base_plot_dir=report_plots_folder,
        plt_name="doublets_umap",
        width=9, height=7
      )
    }
  }

  if (!("Doublet_classifications" %in% names(scrna@meta.data))) {
    message("Doublet_classifications not in meta.data; skipping doublet-specific plots and table.")
    return(invisible(NULL))
  }

  Idents(object = scrna) <- "name"
  feats_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))
  plt1 = VlnPlot(
    object = scrna,
    features = feats_to_plot,
    ncol = 2,
    cols = col_def,
    pt.size = 0
  ) %+safe% ggtitle("All cells")
  plt1_linear = VlnPlot(
    object = scrna,
    features = feats_to_plot,
    ncol = 1,
    cols = col_def,
    pt.size = 0
  ) %+safe% ggtitle("All cells")
  save_ggplot_formats(
    plt=plt1,
    base_plot_dir=report_plots_folder,
    plt_name="vlnplot_doublet_all",
    width=9, height=7
  )

  scrna_subset <- subset(scrna, Doublet_classifications == "Singlet")
  plt2 = VlnPlot(
    object = scrna_subset,
    features = feats_to_plot,
    ncol = 2,
    cols = col_def,
    pt.size=0
  ) %+safe% ggtitle("Singlets")
  plt2_linear = VlnPlot(
    object = scrna_subset,
    features = feats_to_plot,
    ncol = 1,
    cols = col_def,
    pt.size=0
  ) %+safe% ggtitle("Singlets")
  save_ggplot_formats(
    plt=plt2,
    base_plot_dir=report_plots_folder,
    plt_name="vlnplot_doublet_singlets",
    width=9, height=7
  )

  scrna_subset <- subset(scrna, Doublet_classifications == "Doublet")
  plt3 = VlnPlot(
    object = scrna_subset,
    features = feats_to_plot,
    ncol = 2,
    cols = col_def,
    pt.size=0
  )
  plt3_linear = VlnPlot(
    object = scrna_subset,
    features = feats_to_plot,
    ncol = 1,
    cols = col_def,
    pt.size=0
  )
  save_ggplot_formats(
    plt=plt3,
    base_plot_dir=report_plots_folder,
    plt_name="vlnplot_doublet_doublets",
    width=9, height=7
  )
  vln_plt <- patchwork::wrap_plots(plt1_linear, plt2_linear, plt3_linear, ncol = 3)
  save_ggplot_formats(
    plt = vln_plt,
    base_plot_dir=report_plots_folder,
    plt_name = "vlnplot_doublet_combined",
    width=27, height=28
  )

  meta <- scrna@meta.data
  stSample <- meta %>%
    group_by(name) %>%
    summarise(
      Doublets=sum(Doublet_classifications == "Doublet"),
      Singlets=sum(Doublet_classifications == "Singlet")
    )

  save_object(
    stSample,
    file.path(report_tables_folder,"doublets_table.Rds"),
    COMPRESSION_FORMAT
  )
}
