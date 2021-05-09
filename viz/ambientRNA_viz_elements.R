####################################################
# ambient RNA
####################################################
## if load celda here, Seurat methods would be replaced by old version
#suppressPackageStartupMessages(library(celda))

#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(Seurat))


ambientRNA_elements <- function(scrna){
  plt <- list()

  plt[[1]] <- FeaturePlot(scrna, features = "AmbientRNA") + ggtitle(label = "Ambient RNA Contamintaion")

  plt[[2]] <- VlnPlot(object = scrna, features = "AmbientRNA", group.by = "decontX_clusters", pt.size=0.1)

  save_ggplot_formats(
    plt=plt[[1]],
    base_plot_dir=report_plots_folder,
    plt_name="ambient_rna",
    width=9, height=7
  )
  save_ggplot_formats(
    plt=plt[[2]],
    base_plot_dir=report_plots_folder,
    plt_name="ambient_rna_vln",
    width=9, height=7
  )

  meta <- scrna@meta.data
  stSample <- meta %>%
    group_by(name) %>%
    summarise(
      ambientRNA.Mean=mean(AmbientRNA),
      ambientRNA.Median=median(AmbientRNA),
    )
  saveRDS(stSample,file.path(report_tables_folder,"ambientRNA_postfilter.RDS"))
}


