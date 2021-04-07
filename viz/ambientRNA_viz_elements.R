####################################################
# ambient RNA
####################################################
library(celda)
library(ggplot2)
library(dplyr)
library(Seurat)
# savedir <- "/home/grasshoff/labcluster/ambientRNA/results/"
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
}else{
  scrna <- readRDS(file = file.path(savedir, "scrna_phase_preprocess.Rds"))
}
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

meta <- scrna@meta.data
stSample <- meta %>%
  group_by(name) %>%
  summarise(
    ambientRNA.Mean=mean(AmbientRNA),
    ambientRNA.Median=median(AmbientRNA),
  )
saveRDS(stSample,file.path(report_tables_folder,"ambientRNA_postfilter.RDS"))

