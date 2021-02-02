####################################################
# ambient RNA
####################################################
library(celda)
library(ggplot2)
# savedir <- "/home/grasshoff/labcluster/ambientRNA/results/"
scrna <- readRDS(file = file.path(savedir, "scrna_phase_comparing.Rds"))
plt <- list()

plt[[1]] <- FeaturePlot(scrna, features = "soupX_contamination") + ggtitle(label = "Ambient RNA Contamintaion")

plt[[2]] <- VlnPlot(object = scrna, features = "soupX_contamination", group.by = "seurat_clusters")

plt[[3]] <- DimPlot(scrna, group.by = cluster)

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

save_ggplot_formats(
  plt=plt[[3]],
  base_plot_dir=report_plots_folder,
  plt_name="ambient_cluster",
  width=9, height=7
)
