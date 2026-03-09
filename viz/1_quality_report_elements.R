####################################################
# pre filtering
####################################################
quality_report_elements <- function(){

  scrna <- load_object(file_name = file.path(savedir, "scrna_rawdata.Rds"))
  scrna <- safe_join_layers(scrna)  # v5: join layers

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
  # clinical metadata (optional, per-patient barplots)
  ####################################################
  if (!is.null(clinical_meta)) {
    message("### Making clinical metadata barplots")
    meta_samples <- clinical_meta[intersect(names(data_src), rownames(clinical_meta)), , drop=FALSE]
    meta_samples$sample <- factor(rownames(meta_samples), levels = names(data_src))
    meta_samples$stage  <- stage_lst[rownames(meta_samples)]
    sample_col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(meta_samples$sample)))
    names(sample_col_def) <- levels(meta_samples$sample)

    for (col in colnames(clinical_meta)) {
      if (is.numeric(meta_samples[[col]])) {
        plt <- ggplot(meta_samples, aes(x = stage, y = .data[[col]], color = sample, group = sample))

        # Add a clean mean bar per stage behind sample points.
        plt <- plt +
          stat_summary(
            fun = mean,
            fun.min = mean,
            fun.max = mean,
            geom = "crossbar",
            width = 0.62,
            fatten = 0,
            color = "grey25",
            linewidth = 1.15,
            alpha = 0.9,
            na.rm = TRUE
          )

        plt <- plt +
          geom_point(size = 4, position = position_dodge(width = 0.5)) +
          scale_color_manual(values = sample_col_def) +
          theme_minimal() +
          ggtitle(paste("Patient", col)) +
          xlab("") + ylab(col) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        prop_data <- meta_samples %>%
          count(.data[[col]], stage) %>%
          group_by(stage) %>%
          mutate(prop = n / sum(n)) %>%
          ungroup()
        prop_data[[col]] <- as.factor(prop_data[[col]])
        n_stages <- length(unique(prop_data$stage))
        cat_levels <- levels(prop_data[[col]])
        cat_col_def <- ggsci_pal(option = replicates_viridis_opt)(length(cat_levels))
        names(cat_col_def) <- cat_levels
        plt <- ggplot(prop_data, aes(x = "", y = prop, fill = .data[[col]])) +
          geom_bar(stat = "identity", width = 1) +
          geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
          coord_polar("y", start = 0) +
          facet_wrap(~stage, ncol = n_stages) +
          scale_fill_manual(values = cat_col_def) +
          theme_minimal() +
          ggtitle(paste("Proportion of", col, "by stage")) +
          xlab("") + ylab("") +
          labs(fill = col)
      }
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = paste0("qc_clinical_sampleplot_", col),
        width = 9, height = 5
      )
    }

    if (!is.null(extra_stage_cols)) {
      message("### Making clinical metadata dotplots by extra stage groups")
      for (esc_col in names(extra_stage_cols)) {
        esc_def <- extra_stage_cols[[esc_col]]
        esc_meta_col <- if (is.list(esc_def)) esc_col else esc_def
        if (esc_meta_col %in% names(meta_samples)) {
          meta_samples[[esc_meta_col]] <- as.factor(meta_samples[[esc_meta_col]])
          for (clin_col in colnames(clinical_meta)) {
            if (is.numeric(meta_samples[[clin_col]])) {
              agg_data <- meta_samples %>%
                group_by(.data[[esc_meta_col]]) %>%
                summarise(
                  mean_val = mean(.data[[clin_col]], na.rm = TRUE),
                  sd_val = sd(.data[[clin_col]], na.rm = TRUE),
                  n = n(),
                  .groups = 'drop'
                )
              plt <- ggplot(agg_data, aes(x = .data[[esc_meta_col]], y = mean_val, color = .data[[esc_meta_col]])) +
                geom_point(size = 5) +
                geom_errorbar(
                  aes(ymin = mean_val - sd_val / sqrt(n), ymax = mean_val + sd_val / sqrt(n)),
                  width = 0.2
                ) +
                theme_minimal() +
                ggtitle(paste("Mean", clin_col, "by", esc_meta_col)) +
                xlab(esc_meta_col) + ylab(paste("Mean", clin_col))
            } else {
              plt <- ggplot(meta_samples, aes(x = .data[[clin_col]], fill = .data[[esc_meta_col]])) +
                geom_bar(position = "dodge") +
                theme_minimal() +
                ggtitle(paste(clin_col, "by", esc_meta_col)) +
                xlab(clin_col) + labs(fill = esc_meta_col)
            }
            save_ggplot_formats(
              plt = plt,
              base_plot_dir = report_plots_folder,
              plt_name = paste0("clinical_dotplot_", esc_col, "_", clin_col),
              width = 7, height = 5
            )
          }
        }
      }
    }

    save_object(
      meta_samples,
      file.path(report_tables_folder, "clinical_meta_summary.RDS"),
      COMPRESSION_FORMAT
    )
  }

  ####################################################
  # post filtering
  ####################################################
  if(identical(cluster,"singleton")){
    scrna <- load_object(file_name = file.path(savedir, "scrna_phase_singleton.Rds"))
  }else{
    scrna <- load_object(file_name = file.path(savedir, "scrna_phase_preprocess.Rds"))
  }
  #scrna <- safe_join_layers(scrna)  # v5: join layers for downstream GetAssayData / VariableFeatures

  Idents(object = scrna)<- "name"
  if (DOUBLET_SWITCH=='off'){
    feats_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  }else{
    feats_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "pANN")
  }
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
  col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(Idents(scrna))))
  ps <- lapply(feats_to_plot, function(fea){
    VlnPlot(
      object = scrna,
      features = fea,
      group.by="name",
      cols = col_def,
      pt.size=0
    ) %+safe% NoLegend()
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


  if (DOUBLET_SWITCH=='off'){
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
  }else{
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
        pANN.Mean = mean(pANN),
        pANN.Median = median(pANN),
        Cells = sum(cells)
      )
  }


  save_object(
    stSample,
    file.path(report_tables_folder,"stSample_postfilter.RDS"),
    COMPRESSION_FORMAT
  )

  if (DOUBLET_SWITCH=='off'){
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
  }else{
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
        pANN.Mean = mean(pANN),
        pANN.Median = median(pANN),
        Cells = sum(cells)
      )
  }


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

  if (DOUBLET_SWITCH =='off'){
    plt = patchwork::wrap_plots(list(p1, p2, p3), ncol=1)
  }else{
    p4 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "pANN", cols=col_def)
    plt = patchwork::wrap_plots(list(p1, p2, p3, p4), ncol=1)
  }


  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name="featurescatter_ncountRNA",
    width=9, height=18
  )


  ## High variable genes (Seurat v5: ensure HVF metadata exists before VariableFeaturePlot)
  col_def <- c(base_color, pos_color)
  DefaultAssay(scrna) <- "RNA"
  vf <- VariableFeatures(scrna)
  if (length(vf) == 0) {
    scrna <- NormalizeData(scrna, verbose = FALSE)
    scrna <- FindVariableFeatures(scrna, verbose = FALSE)
    vf <- VariableFeatures(scrna)
  }
  if (length(vf) > 0) {
    top10 <- head(vf, 10)
    plot1 <- VariableFeaturePlot(scrna, cols = col_def)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plt <- patchwork::wrap_plots(list(plot1, plot2))
    save_ggplot_formats(
      plt = plt,
      base_plot_dir = report_plots_folder,
      plt_name = "high_var_genes",
      width = 13, height = 5
    )
  }

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
