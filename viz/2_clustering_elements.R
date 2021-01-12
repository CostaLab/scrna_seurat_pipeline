## Clusters Resolution
if(identical(cluster,"singleton")){
  # scrna <- readRDS(file=file.path(savedir, "scrna_singleton_markergenes.Rds"))
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_singleton.Rds"))
}else{
  scrna <- readRDS(file=file.path(savedir, "scrna_phase_comparing.Rds"))
}

for(cluster_use in available_clusters){

  message(paste0("### Producing elements for cluster: ",cluster_use))

  fisher_cluster_name <- paste0("fishertest_", cluster_use)

  if(cluster_use %ni% names(scrna@meta.data)){
    stop(glue("ERROR:There's no this {cluster_use} slot, please check!!!"))
  }

  if(fisher_cluster_name %ni% names(scrna@tools)){
    stop(glue("ERROR:fishertest hasn't been calculated for cluster {cluster_use}\n Please run [scrna_fishertest_clusters]!!!"))
  }

  pref_def = "integrated_snn_res."
  if(cluster_use == "harmony_inte_clusters") pref_def = "RNA_snn_res."
  if(cluster_use == "singleton") pref_def = "RNA_snn_res."

  umap_reduction = "DEFAULT_UMAP"
  if(cluster_use == "harmony_inte_clusters") umap_reduction = "harmony_UMAP"
  if(cluster_use == "seurat_inte_clusters") umap_reduction = "INTE_UMAP"
  if(cluster_use == "singleton") umap_reduction = "SINGLE_UMAP"

  message("### Making cluster tree")
  plt = clustree(scrna, prefix = pref_def)
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("clustree_resolution_",cluster_use),
    width=13, height=10
  )


  message("### Making umap resolution list")
  # UMAP resolution list
  cluster_de_list <- scrna@tools$de_batch
  # FIXME should be the resolution vector the user defined right?
  names(cluster_de_list) = as.character(seq(0.1, 0.8, 0.1))
  nms <- names(cluster_de_list)
  plist <- list()
  for(nm in nms){

    plist[[nm]] <- DimPlot(
      scrna,
      reduction = umap_reduction,
      group.by =  paste0(pref_def, nm), #FIXME group.by remains the same regardless of reduction?
      label=T, label.size=8
    ) + ggtitle(sprintf("resolution %s", nm))

  }

  plt = patchwork::wrap_plots(plist, ncol=2)
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_resolution_list_",cluster_use),
    width=15, height=20
  )

  message("### Making umap grouped by name")
  ## Clusters
  plt = DimPlot(scrna, reduction = umap_reduction, group.by = "name", cols=colours)
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_groupby_name_",cluster_use),
    width=9, height=7
  )
  message("### Making umap grouped by clusters")
  plt = DimPlot(scrna, reduction = umap_reduction, group.by = cluster_use, label=T, label.size=8)
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_groupby_cluster_",cluster_use),
    width=9, height=7
  )

  saveRDS(
    table(scrna$name, scrna@meta.data[, cluster_use]),
    file.path(report_tables_folder, paste0("stSample_table_cluster_",cluster_use,".RDS"))
  )
  saveRDS(
    table(scrna$stage, scrna@meta.data[, cluster_use]),
    file.path(report_tables_folder, paste0("stCond_table_cluster_",cluster_use,".RDS"))
  )

  if(cluster_use != "singleton"){
    ## Clusters Statistics
    message("### Making cluster statistics barplot")
    df <- scrna@tools[[fisher_cluster_name]]

    if(all.is.numeric(df$Cluster)){ ## set int order if all cluster name are integers
      df$Cluster <- factor(
        as.character(df$Cluster),
        levels = sort(as.numeric(as.character(df$Cluster)), decreasing=T)
      )
    }

    shift_trans = function(d = 0) {
      scales::trans_new(
        "shift",
        transform = function(x) x - d,
        inverse = function(x) x + d)
    }

    if(length(table(scrna$stage))<=2){
      nm <- scrna@tools[["meta_order"]][["stage"]][1]
      nm2 <- scrna@tools[["meta_order"]][["stage"]][2]

      aodds <- glue("odds.ratio_{nm}.vs.{nm2}", )
      plt <- ggplot(data=df, aes_string(x = "Cluster", y = aodds, fill = "Cluster")) +
        geom_bar(stat="identity") +
        coord_flip() +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_y_continuous(trans = shift_trans(1)) + geom_text(data = df,
                    aes_string("Cluster", 1, label = glue("pval.adjust_{nm}.vs.{nm2}")),
                    position = "identity",
                    size=4) +
        ggtitle(sprintf("%s vs %s (odds ratio > 1 means more %s)", nm, nm2, nm))+
        xlab("Cluster") +
        ylab("odds ratio")+
        theme_minimal()

      save_ggplot_formats(
        plt=plt,
        base_plot_dir=report_plots_folder,
        plt_name=paste0("cluster_statistics_",nm,"-vs-",nm2,"_",fisher_cluster_name),
        width=9, height=7
      )

    }else{
      pairs <- combn(1:length(table(scrna$stage)), 2)
      n <- length(pairs)/2
      for (i in 1:n){
        i1 <- pairs[1:2, i][1]
        i2 <- pairs[1:2, i][2]
        nm <- scrna@tools[["meta_order"]][["stage"]][i1]
        nm2 <- scrna@tools[["meta_order"]][["stage"]][i2]

        aodds <- glue("odds.ratio_{nm}.vs.{nm2}")
        plt <- ggplot(data=df, aes_string(x = "Cluster", y = aodds, fill = "Cluster")) +
        geom_bar(stat="identity") +
        coord_flip() +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_y_continuous(trans = shift_trans(1)) +
        geom_text(data = df,
                  aes_string("Cluster", 1, label = glue("pval.adjust_{nm}.vs.{nm2}")),
                  position = "identity",
                  size=4) +
        ggtitle(glue("{nm} vs {nm2} (odds ratio > 1 means more {nm})"))+
        xlab("Cluster") +
        ylab("odds ratio")+
        theme_minimal()
        save_ggplot_formats(
          plt=plt,
          base_plot_dir=report_plots_folder,
          plt_name=paste0("cluster_statistics_",nm,"-vs-",nm2,"_",fisher_cluster_name),
          width=9, height=7
        )
      }
    }


    ## Proportion
    message("### Making sample proportions per cluster barplot")
    name_len <- length(table(scrna@meta.data$name))
    help_sort_func <- ifelse(
      all.is.numeric(unique(scrna@meta.data[, cluster_use])),
      function(x) as.numeric(x) - 1,
      as.character
    )
    scrna@meta.data[,cluster_use] <- help_sort_func(scrna@meta.data[,cluster_use])

    cluster_propotion_before = t(prop.table(x = table(scrna@meta.data$name, scrna@meta.data[, cluster_use]), margin = 2))
    cluster_propotion_before = cluster_propotion_before[,1:name_len]
    cluster_propotion_before_sort = cluster_propotion_before[order(cluster_propotion_before[,1],
                                            cluster_propotion_before[,2],decreasing=TRUE),]

    # barplot uses default graphic device hence we use the std functions to save it
    png(
      filename=file.path(report_plots_folder_png,paste0("cluster_sample_proportions_",cluster_use,".png")),
      width=13,
      height=7,
      units="in",
      type="cairo",
      res=300
    )
    par(xpd=TRUE)
    par(mar = c(5,2,2,12))
    barplot(
      t(cluster_propotion_before_sort), col=colours, xlab="cluster",
      legend.text = colnames(cluster_propotion_before_sort),
      args.legend = list(x ='right', bty='n', inset=c(-0.13,0), xpd = TRUE),
      names.arg = rownames(cluster_propotion_before_sort),
      main = "Propotion of dataset"
    )
    dev.off()

    pdf(
      file=file.path(report_plots_folder_pdf,paste0("cluster_sample_proportions_",cluster_use,".pdf")),
      width=13,
      height=7
    )
    par(xpd=TRUE)
    par(mar = c(5,2,2,12))
    barplot(
      t(cluster_propotion_before_sort), col=colours, xlab="cluster",
      legend.text = colnames(cluster_propotion_before_sort),
      args.legend = list(x ='right', bty='n', inset=c(-0.13,0), xpd = TRUE),
      names.arg = rownames(cluster_propotion_before_sort),
      main = "Propotion of dataset"
    )
    dev.off()



    ## Amount Distribution
    message("### Making sample amount distribution barplot")
    #myel <- SetAllIdent(scrna, "name")

    plt = BarPlot(
      scrna@meta.data[, cluster_use], fill = scrna$name, xlab = "Cluster",
      legend.title = "Replicate", main = "Amount of samples in each cluster")
    plt = plt + scale_fill_manual(values=colours)
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("cluster_sample_distribution_",cluster_use),
      width=14, height=9
    )

  }

  ## Cell cycle phase
  message("### Making cellcycle phases umap")
  plt = DimPlot(scrna, reduction = umap_reduction, group.by = "Phase")
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_ccycle_phase_",cluster_use),
    width=9, height=7
  )

  ## FeaturePlot
  message("### Making umap featureplot with QC elements")
  plt = FeaturePlot(
    scrna,
    features = c("percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA"),
    cols = c("lightgrey", "red"),
    reduction=umap_reduction,
    ncol = 2
  )
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_featureplot_qc_",cluster_use),
    width=9, height=7
  )
  message("### Making umap featureplot with ccycle elements")
  plt = FeaturePlot(
    scrna,
    features = c("CC.Difference","G1.Score", "S.Score", "G2M.Score"),
    cols = c("lightgrey", "red"),
    reduction=umap_reduction,
    ncol = 2
  )
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("umap_featureplot_ccycle_",cluster_use),
    width=9, height=7
  )


  ## Violin Plot
  message("### Making violinplot with QC and ccycle elements")
  plt = VlnPlot(
    scrna,
    features = c(
      "percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA",
      "CC.Difference", "G1.Score", "S.Score", "G2M.Score"
    ),
    group.by = cluster_use,
    pt.size = 0
  )
  save_ggplot_formats(
    plt=plt,
    base_plot_dir=report_plots_folder,
    plt_name=paste0("vlnplot_qc_ccycle_groupby-cluster_",cluster_use),
    width=9, height=7
  )



  ## MCA annotation
  if("MCA_annotate" %in% names(scrna@meta.data)){
    message("### Making umap with MCA annotation")
    plt <- DimPlot(
      object = scrna,
      reduction = umap_reduction,
      pt.size = 0.2,
      group.by = "MCA_annotate"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(colour="blue", size=4, face="bold"),
      legend.text = element_text(size = 7)
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("mca_annotate_",cluster_use),
      width=9, height=7
    )

    # also save plot and information as rds so that it can be later rendered in the report as plotly
    message("### Saving MCA annotation data to produce plotly in report")
    saveRDS(plt,file.path(report_plots_folder,paste0("mca_annotate_plt_",cluster_use,".RDS")))
    saveRDS(
      FetchData(object = scrna, vars = c("MCA_annotate", cluster_use)),
      file.path(report_plots_folder,paste0("mca_annotate_info_",cluster_use,".RDS"))
    )
  }

  ## HCL annotation
  if("HCL_annotate" %in% names(scrna@meta.data)){
    message("### Making umap with HCL annotation")
    plt <- DimPlot(
      object = scrna,
      reduction = umap_reduction,
      pt.size = 0.2,
      group.by = "HCL_annotate"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(colour="blue", size=4, face="bold"),
      legend.text = element_text(size = 7)
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("hcl_annotate_",cluster_use),
      width=9, height=7
    )

    # also save plot and information as rds so that it can be later rendered in the report as plotly
    message("### Saving HCL annotation data to produce plotly in report")
    saveRDS(plt,file.path(report_plots_folder,paste0("hcl_annotate_plt_",cluster_use,".RDS")))
    saveRDS(
      FetchData(object = scrna, vars = c("HCL_annotate", cluster_use)),
      file.path(report_plots_folder,paste0("hcl_annotate_info_",cluster_use,".RDS"))
    )
  }

  ## External Annotation
  if("external_annotation" %in% names(scrna@meta.data)){
    message("### Making umap with external annotation")
    plt <- DimPlot(
      object = scrna,
      reduction = umap_reduction,
      pt.size = 0.2,
      group.by = "external_annotation"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(colour="blue", size=4, face="bold"),
      legend.text = element_text(size = 7)
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("external_annotation_",cluster_use),
      width=9, height=7
    )
    message("### Saving external annotation data to produce plotly in report")
    saveRDS(plt,file.path(report_plots_folder,paste0("external_annotation_plt_",cluster_use,".RDS")))
    saveRDS(
      FetchData(object = scrna, vars = c("external_annotation", cluster_use)),
      file.path(report_plots_folder,paste0("external_annotation_info_",cluster_use,".RDS"))
    )
  }
}
