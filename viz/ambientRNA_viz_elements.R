####################################################
# ambient RNA
####################################################
if(identical(cluster,"singleton")){
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
}else{
  scrna <- readRDS(file = file.path(savedir, "scrna_phase_clustering.Rds"))
}
plt <- list()

plt[[1]] <- FeaturePlot(scrna, features = "AmbientRNA") + ggtitle(label = "Ambient RNA Contamintaion")

plt[[2]] <- VlnPlot(object = scrna, features = "AmbientRNA", group.by = "decontX_clusters", pt.size=0.1)

plt[[3]] <- VlnPlot(object = scrna, features = "AmbientRNA", group.by = "decontX_clusters")

plt[[4]] <- VlnPlot(object = scrna, features = "AmbientRNA_Harmony", group.by = "decontX_clusters_harmony")

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
  plt_name="ambient_rna_vln",
  width=9, height=7
)
save_ggplot_formats(
  plt=plt[[4]],
  base_plot_dir=report_plots_folder,
  plt_name="ambient_rna_vln_harmony",
  width=9, height=7
)

meta <- scrna@meta.data
stSample <- meta %>%
  group_by(name) %>%
  summarise(
    ambientRNA.Mean=mean(AmbientRNA),
    ambientRNA.Median=median(AmbientRNA),
    ambientRNA.Harmony.Mean=mean(AmbientRNA_Harmony),
    ambientRNA.Harmony.Median=median(AmbientRNA_Harmony)
  )
saveRDS(stSample,file.path(report_tables_folder,"ambientRNA_postfilter.RDS"))
