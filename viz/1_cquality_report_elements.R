####################################################
# QC on the current clusters
####################################################
cquality_report_elements <- function(scrna){

  Idents(object = scrna)<- "name"

  feats_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))

  plt = VlnPlot(
    object = scrna,
    features = feats_to_plot,
    ncol=2,
    cols = col_def,
    pt.size=0
  )
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="postfilter_vlnplot",
    width=9, height=14
  )

  feats_to_plot = c("G1.Score", "S.Score", "G2M.Score")
  if ("G1.Score" %ni% names(scrna@meta.data)){
    scrna$G1.Score = 1 - scrna$S.Score - scrna$G2M.Score
  }
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))
  ps <- lapply(feats_to_plot, function(fea){
    VlnPlot(
      object = scrna,
      features = fea,
      group.by="name",
      cols = col_def,
      pt.size=0
    ) + NoLegend()
  })

  plt = plot_grid(plotlist=ps, ncol=2)

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="postfilter_cellphase",
    width=9, height=7
  )


  meta <- scrna@meta.data
  meta$cells <- 1

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
      pctRb.Median = median(percent.ribo),
      Cells = sum(cells)
    )
  save_object(
    stSample,
    file.path(report_tables_folder,"stSample_postfilter.RDS"),
    COMPRESSION_FORMAT
  )

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
      pctRb.Median = median(percent.ribo),
      Cells = sum(cells)
    )
  save_object(
    stCond,
    file.path(report_tables_folder,"stCond_postfilter.RDS"),
    COMPRESSION_FORMAT
  )



  # qc feature scatterplots
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))
  p1 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "percent.mt", cols=col_def)
  p2 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "percent.ribo", cols=col_def)
  p3 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=col_def)

  plt = patchwork::wrap_plots(list(p1, p2, p3), ncol=1)

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="featurescatter_ncountRNA",
    width=9, height=13
  )

}
