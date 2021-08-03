#doublet_switch###################################################
# doublets 
####################################################
doubletdetection_viz_elements <- function(scrna){
  # We have three options.
  # If doublet_switch == "display", we have all the cells still in the Seurat object.
  # We can plot the doublets in the INTE_UMAP.
  # If doublet_switch == "on", we have removed the doublets from the Seurat object.
  # We can use the output from the doublet detection function (scrna_DoubletAnnotated.Rds) and plot the UMAP we created for this purpose.
  # If doublet_switch == "off", we just put out a message stating that doublet detection was not performed.
  if(doublet_switch == "display"){
    scrna <- readRDS(file.path(savedir, "scrna_phase_comparing.Rds"))
    plt <- DimPlot(scrna, group.by = "Doublet_classifications")
    save_ggplot_formats(
      plt=plt[[1]],
      base_plot_dir=report_plots_folder,
      plt_name="doublets",
      width=9, height=7
    )
    meta <- scrna@meta.data
    stSample <- meta %>%
      group_by(name) %>%
      summarise(
        Doublets=sum(Doublet_classifications == "Doublet"),
        Singlets=sum(Doublet_classifications == "Singlet")
      )
    saveRDS(stSample,file.path(report_tables_folder,"doublets_table.Rds"))
  }
  if(doublet_switch == "on"){
    scrna <- readRDS(file.path(savedir, "scrna_DoubletAnnotated.Rds"))
    plt	<- DimPlot(scrna, group.by = "Doublet_classifications", reduction = "DOUBLET_UMAP")
    save_ggplot_formats(
      plt=plt[[1]],
      base_plot_dir=report_plots_folder,
      plt_name="doublets",
      width=9, height=7
    )
    meta <- scrna@meta.data
    stSample <- meta %>%
      group_by(name) %>%
      summarise(
        Doublets=sum(Doublet_classifications == "Doublet"),
        Singlets=sum(Doublet_classifications == "Singlet")
    )
    saveRDS(stSample,file.path(report_tables_folder,"doublets_table.Rds"))
  }
}
