####################################################
# pre filtering
####################################################
quality_report_elements <- function(){

  scrna <- load_object(file_name = file.path(savedir, "scrna_rawdata.Rds"))

  Idents(object = scrna) <- "name"

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
    plt_name="prefilter_vlnplot",
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
    file.path(report_tables_folder,"stSample_prefilter.RDS"),
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
    file.path(report_tables_folder,"stCond_prefilter.RDS"),
    COMPRESSION_FORMAT
  )

  ####################################################
  # post filtering
  ####################################################
  if(identical(cluster,"singleton")){
    scrna <- load_object(file_name = file.path(savedir, "scrna_phase_singleton.Rds"))
  }else{
    scrna <- load_object(file_name = file.path(savedir, "scrna_phase_preprocess.Rds"))
  }

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
    width=9, height=7
  )


  feats_to_plot = c("G1.Score", "S.Score", "G2M.Score")
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
    width=9, height=18
  )


  ## High variable genes
  col_def <- c(base_color,pos_color)
  top10 <- head(VariableFeatures(scrna), 10)
  plot1 <- VariableFeaturePlot(scrna,cols = col_def)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

  plt = patchwork::wrap_plots(list(plot1, plot2))

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="high_var_genes",
    width=13, height=5
  )

  ## Cellcycle scaling
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))

  ### before
  plt = DimPlot(
    scrna,
    reduction="BCELLCYCLE_PCA",
    cols=col_def
  )

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="cellcycle_scaling_before",
    width=10, height=8
  )

  ### after
  plt = DimPlot(
    scrna,
    reduction="CELLCYCLED_PCA",
    cols=col_def
  )

  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="cellcycle_scaling_after",
    width=10, height=8
  )
}

