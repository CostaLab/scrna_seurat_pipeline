####################################################
# ambient RNA soupx
####################################################
if(identical(cluster,"singleton")){
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
}else{
  scrna <- readRDS(file = file.path(savedir, "scrna_phase_clustering.Rds"))
}
meta <- scrna@meta.data
stSample <- meta %>%
  group_by(name) %>%
  summarise(
    ambientRNA.Mean=mean(soupx_contamination),
    ambientRNA.Median=median(soupx_contamination),
  )
saveRDS(stSample,file.path(report_tables_folder,"ambientRNA_soupx_postfilter.Rds"))


plt <- list()

plt[[1]] <- FeaturePlot(scrna, features = "soupx_contamination") + ggtitle(label = "Ambient RNA Contamintaion SoupX")

plt[[2]] <- VlnPlot(object = scrna, features = "soupx_contamination", group.by = cluster, pt.size=0.1)

plt[[3]] <- VlnPlot(object = scrna, features = "soupx_contamination", group.by = "name", pt.size=0.1)

save_ggplot_formats(
  plt=plt[[1]],
  base_plot_dir=report_plots_folder,
  plt_name="ambient_rna_soupx",
  width=9, height=7
)
save_ggplot_formats(
  plt=plt[[2]],
  base_plot_dir=report_plots_folder,
  plt_name="ambient_rna_soupx_vln",
  width=9, height=7
)
save_ggplot_formats(
  plt=plt[[3]],
  base_plot_dir=report_plots_folder,
  plt_name="ambient_rna_soupx_sample_vln",
  width=9, height=7
)

# We generate the barplots.
# We read the sheets of the xlsx files. For each sample there should be a file.
samples <- unique(scrna$name)
background <- list()
for(i in 1:length(samples)){
  background[[i]] <- read.xlsx(file.path(report_data_folder, "BackgroundEstimation.xlsx"), sheet = samples[i])
  background[[i]] <- head(background[[i]][order(background[[i]]$est, decreasing = TRUE),], n = 20)
  background[[i]]$gene <- factor(background[[i]]$gene, levels = background[[i]]$gene)
  p <- ggplot(background[[i]], aes(x = gene,y = est)) + geom_bar(stat="identity", fill = "steelblue") +
       ylab("fraction of soup reads") +
       theme_minimal() + ggtitle(paste0("SoupX Background Estimation ", samples[i])) +
       theme(axis.text.x = element_text(angle = 45))

  save_ggplot_formats(
    plt=p,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("soupx_background_", samples[i]),
    width=9, height=7
  )
}
