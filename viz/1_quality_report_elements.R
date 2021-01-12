####################################################
# pre filtering
####################################################
scrna <- readRDS(file = file.path(savedir, "scrna_rawdata.Rds"))

Idents(object = scrna) <- "name"

plt = VlnPlot(
  object = scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol=2,
  cols = colours,
  pt.size=0
)
save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="prefilter_vlnplot",
  width=9, height=7
)

meta <- scrna@meta.data

stSample <- meta %>%
  group_by(name) %>%
  summarise(
    nCount.Mean=mean(nCount_RNA),
    nCount.Median=median(nCount_RNA),
    nFeature.Mean=mean(nFeature_RNA),
    nFeature.Median=median(nFeature_RNA),
    pctMt.Mean=mean(percent.mt),
    pctMt.Median=median(percent.mt),
    pctRb.Mean = mean(percent.ribo),
    pctRb.Median = median(percent.ribo)
  )
saveRDS(stSample,file.path(report_tables_folder,"stSample_prefilter.RDS"))

stCond <- meta %>%
  group_by(stage) %>%
  summarise(
    nCount.Mean=mean(nCount_RNA),
    nCount.Median=median(nCount_RNA),
    nFeature.Mean=mean(nFeature_RNA),
    nFeature.Median=median(nFeature_RNA),
    pctMt.Mean=mean(percent.mt),
    pctMt.Median=median(percent.mt),
    pctRb.Mean = mean(percent.ribo),
    pctRb.Median = median(percent.ribo)
  )
saveRDS(stCond,file.path(report_tables_folder,"stCond_prefilter.RDS"))


####################################################
# post filtering
####################################################
if(identical(cluster,"singleton")){
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
}else{
  scrna <- readRDS(file = file.path(savedir, "scrna_phase_preprocess.Rds"))
}

Idents(object = scrna)<- "name"

plt = VlnPlot(
  object = scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol=2,
  cols = colours,
  pt.size=0
)
save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="postfilter_vlnplot",
  width=9, height=7
)


features = c("G1.Score", "S.Score", "G2M.Score")
ps <- lapply(features, function(fea){
             VlnPlot(object = scrna,
                  features = fea,
                  group.by="name",
                  cols = colours,
                  pt.size=0) + NoLegend()

})

plt = plot_grid(plotlist=ps, ncol=2)

save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="postfilter_cellphase",
  width=9, height=7
)


meta <- scrna@meta.data

stSample <- meta %>%
  group_by(name) %>%
  summarise(
    nCount.Mean=mean(nCount_RNA),
    nCount.Median=median(nCount_RNA),
    nFeature.Mean=mean(nFeature_RNA),
    nFeature.Median=median(nFeature_RNA),
    pctMt.Mean=mean(percent.mt),
    pctMt.Median=median(percent.mt),
    pctRb.Mean = mean(percent.ribo),
    pctRb.Median = median(percent.ribo)
  )
saveRDS(stSample,file.path(report_tables_folder,"stSample_postfilter.RDS"))

stCond <- meta %>%
  group_by(stage) %>%
  summarise(
    nCount.Mean=mean(nCount_RNA),
    nCount.Median=median(nCount_RNA),
    nFeature.Mean=mean(nFeature_RNA),
    nFeature.Median=median(nFeature_RNA),
    pctMt.Mean=mean(percent.mt),
    pctMt.Median=median(percent.mt),
    pctRb.Mean = mean(percent.ribo),
    pctRb.Median = median(percent.ribo)
  )
saveRDS(stCond,file.path(report_tables_folder,"stCond_postfilter.RDS"))



# qc feature scatterplots
p1 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 ="percent.mt", cols=colours)
p2 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 ="percent.ribo",cols=colours)
p3 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols=colours)

plt = patchwork::wrap_plots(list(p1, p2, p3), ncol=1)

save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="featurescatter_ncountRNA",
  width=9, height=18
)


## High variable genes
top10 <- head(VariableFeatures(scrna), 10)
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plt = patchwork::wrap_plots(list(plot1, plot2))

save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="high_var_genes",
  width=13, height=5
)

## Cellcycle scaling
### before
plt = DimPlot(scrna, reduction="BCELLCYCLE_PCA")

save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="cellcycle_scaling_before",
  width=10, height=8
)

### after
plt = DimPlot(scrna, reduction="CELLCYCLED_PCA")

save_ggplot_formats(
  plt=plt,
  base_plot_dir=report_plots_folder,
  plt_name="cellcycle_scaling_after",
  width=10, height=8
)
