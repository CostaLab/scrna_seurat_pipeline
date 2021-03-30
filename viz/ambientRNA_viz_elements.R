####################################################
# ambient RNA
####################################################
library(celda)
library(ggplot2)
# savedir <- "/home/grasshoff/labcluster/ambientRNA/results/"
scrna <- readRDS(file = file.path(savedir, "scrna_phase_preprocess.Rds"))
plt <- list()

plt[[1]] <- FeaturePlot(scrna, features = "AmbientRNA") + ggtitle(label = "Ambient RNA Contamintaion")

plt[[2]] <- VlnPlot(object = scrna, features = "AmbientRNA", group.by = "decontX_clusters")

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


