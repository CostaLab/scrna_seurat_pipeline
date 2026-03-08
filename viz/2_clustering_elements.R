## Clusters Resolution

clustering_elements <- function(scrna){
  for(cluster_use in available_clusters){

    message(paste0("### Producing elements for cluster: ",cluster_use))

    fisher_cluster_name <- paste0("fishertest_", cluster_use)

    if(cluster_use %ni% names(scrna@meta.data)){
      stop(glue("ERROR:There's no this {cluster_use} slot, please check!!!"))
    }

    if(fisher_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:fishertest hasn't been calculated for cluster {cluster_use}\n Please run [scrna_fishertest_clusters]!!!"))
    }


    umap_reduction = "DEFAULT_UMAP"
    if(cluster_use == "harmony_inte_clusters") umap_reduction = "harmony_UMAP"
    if(cluster_use == "seurat_inte_clusters") umap_reduction = "INTE_UMAP"
    if(cluster_use == "singleton") umap_reduction = "SINGLE_UMAP"


    message("### Making umap grouped by name")
    group_by <- "name"
    col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(scrna@meta.data[,group_by])))
    ## Clusters
    plt = DimPlot(
      scrna,
      reduction = umap_reduction,
      group.by = group_by,
      cols=col_def
    )

    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_groupby_name_",cluster_use),
      width=9, height=7
    )
    message("### Making umap grouped by clusters")
    group_by <- cluster_use
    col_def <- rev(ggsci_pal(option = cluster_viridis_opt)(length(unique(scrna@meta.data[,group_by]))))
    plt = DimPlot(
      scrna,
      reduction = umap_reduction,
      group.by = group_by,
      cols=col_def,
      label=TRUE,
      label.size=8
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_groupby_cluster_",cluster_use),
      width=9, height=7
    )

    message("### Making umap grouped by stage")
    group_by <- "stage"
    col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(scrna@meta.data[,group_by])))
    plt = DimPlot(
      scrna,
      reduction = umap_reduction,
      group.by = group_by,
      cols=col_def
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_groupby_stage_",cluster_use),
      width=9, height=7
    )

    ## Clinical metadata: violin plots & UMAPs per cluster
    if (!is.null(clinical_meta)) {
      numeric_cols <- colnames(clinical_meta)[sapply(clinical_meta, is.numeric)]
      cat_cols <- colnames(clinical_meta)[sapply(clinical_meta, function(x) is.character(x) || is.factor(x))]

      for (col in numeric_cols) {
        if (col %in% colnames(scrna@meta.data)) {
          message(paste0("### Making clinical barplot: ", col))
          meta_cluster <- scrna@meta.data[, c(cluster_use, col)]
          meta_cluster[[col]] <- as.factor(meta_cluster[[col]])
          # Reorder factor levels by numeric value for proper stacking order
          if (all(!is.na(suppressWarnings(as.numeric(levels(meta_cluster[[col]])))))) {
            meta_cluster[[col]] <- factor(meta_cluster[[col]], levels = sort(as.numeric(levels(meta_cluster[[col]]))))
          }
          prop_data <- meta_cluster %>%
            count(.data[[cluster_use]], .data[[col]]) %>%
            group_by(.data[[cluster_use]]) %>%
            mutate(prop = n / sum(n)) %>%
            ungroup()
          # Check if the column values are numeric (for color gradient)
          is_numeric_col <- all(!is.na(suppressWarnings(as.numeric(levels(meta_cluster[[col]])))))
          n_lvl <- length(unique(meta_cluster[[col]]))
          col_def <- ggsci_pal(option = replicates_viridis_opt)(max(n_lvl, 2))
          if (is_numeric_col) {
            num_vals <- sort(as.numeric(levels(meta_cluster[[col]])))
            names(col_def) <- as.character(num_vals)
            prop_data[[col]] <- factor(prop_data[[col]], levels = as.character(num_vals))
            plt <- ggplot(prop_data, aes(x = .data[[cluster_use]], y = prop, fill = .data[[col]])) + geom_col(position = "fill") + scale_fill_manual(values = col_def, breaks = as.character(num_vals)) + theme_minimal() + ggtitle(paste(col, "proportion by cluster")) + xlab(cluster_use) + ylab("Proportion") + labs(fill = col)
          } else {
            plt <- ggplot(prop_data, aes(x = .data[[cluster_use]], y = prop, fill = .data[[col]])) +
              geom_col(position = "fill") +
              theme_minimal() +
              ggtitle(paste(col, "proportion by cluster")) +
              xlab(cluster_use) + ylab("Proportion") +
              labs(fill = col)
          }
          save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder,
            plt_name = paste0("clustering_clinical_clusterprop_", col, "_", cluster_use),
            width = 10, height = 6)

          message(paste0("### Making clinical UMAP: ", col))
          plt <- FeaturePlot(scrna, features = col, reduction = umap_reduction)
          save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder,
            plt_name = paste0("clinical_umap_", col, "_", cluster_use),
            width = 9, height = 7)
        }
      }

      for (col in cat_cols) {
        if (col %in% colnames(scrna@meta.data)) {
          message(paste0("### Making clinical UMAP (categorical): ", col))
          n_lvl <- length(unique(scrna@meta.data[, col]))
          col_def_cat <- ggsci_pal(option = replicates_viridis_opt)(max(n_lvl, 2))
          plt <- DimPlot(scrna, group.by = col, reduction = umap_reduction, cols = col_def_cat)
          save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder,
            plt_name = paste0("clinical_umap_", col, "_", cluster_use),
            width = 9, height = 7)
        }
      }
    }

    ## Extra stage columns: UMAP + composition barplots
    esc <- scrna@tools[["extra_stage_cols"]]
    if (!is.null(esc) && length(esc) > 0) {
      for (col_name in names(esc)) {
        defn <- esc[[col_name]]
        meta_col <- if (is.list(defn)) col_name else defn
        if (!(meta_col %in% colnames(scrna@meta.data))) next
        groups <- sort(na.omit(unique(as.character(scrna@meta.data[, meta_col]))))
        if (length(groups) < 2) next

        message(paste0("### Making extra stage UMAP: ", col_name))
        n_lvl <- length(groups)
        esc_cols <- ggsci_pal(option = replicates_viridis_opt)(max(n_lvl, 2))
        plt <- DimPlot(scrna, group.by = meta_col, reduction = umap_reduction, cols = esc_cols)
        save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder,
          plt_name = paste0("extrastage_umap_", col_name, "_", cluster_use),
          width = 9, height = 7)

        message(paste0("### Making extra stage composition: ", col_name))
        comp_tbl <- table(scrna@meta.data[, meta_col], scrna@meta.data[, cluster_use])
        comp_df <- as.data.frame(prop.table(comp_tbl, margin = 2))
        colnames(comp_df) <- c("Group", "Cluster", "Proportion")
        plt <- ggplot(comp_df, aes(x = Cluster, y = Proportion, fill = Group)) +
          geom_col(position = "fill") +
          theme_minimal() +
          ggtitle(paste("Cluster composition by", col_name)) +
          xlab("Cluster") + ylab("Proportion") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder,
          plt_name = paste0("extrastage_composition_", col_name, "_", cluster_use),
          width = 10, height = 6)

        fisher_key <- sprintf("fisher_extrastage_%s_%s", col_name, cluster_use)
        fisher_data <- scrna@tools[[fisher_key]]
        if (!is.null(fisher_data) && !is.null(fisher_data$fisher)) {
          message(paste0("### Making extra stage Fisher heatmap: ", col_name))
          or_vals <- sapply(fisher_data$fisher, function(x) if (!is.null(x$estimate)) x$estimate else NA)
          pvals <- sapply(fisher_data$fisher, function(x) x$p.value)
          fisher_df <- data.frame(Cluster = names(or_vals), OddsRatio = or_vals, pvalue = pvals, stringsAsFactors = FALSE)
          fisher_df$log2_OR <- log2(fisher_df$OddsRatio)
          fisher_df$pval_label <- ifelse(fisher_df$pvalue < 0.001, sprintf("%.1e", fisher_df$pvalue), sprintf("%.3f", fisher_df$pvalue))
          fisher_col_def <- rev(ggsci_pal(option = cluster_viridis_opt)(length(unique(fisher_df$Cluster))))
          # Extract comparison groups from col_name (e.g., "sex_M.vs.F" -> positive = M)
          pos_group <- gsub(".*_([^_]+).vs.*", "\\1", col_name)
          plt <- ggplot(fisher_df, aes(x = Cluster, y = log2_OR, fill = factor(Cluster, levels = sort(as.numeric(Cluster))))) + geom_col() + geom_hline(yintercept = 0, linetype = "dashed", color = "black") + scale_fill_manual(values = fisher_col_def) + geom_text(aes(y = 0, label = pval_label), vjust = ifelse(fisher_df$log2_OR >= 0, -0.5, 1.5), color = "black", size = 3) + theme_minimal() + ggtitle(paste0("Fisher test: ", col_name, " vs cluster (log2 Odds Ratio, positive=", pos_group, ")")) + xlab("Cluster") + ylab("log2(Odds Ratio)") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill = "Cluster")
          save_ggplot_formats(plt = plt, base_plot_dir = report_plots_folder, plt_name = paste0("extrastage_fisher_", col_name, "_", cluster_use), width = 10, height = 6)
        }

        save_object(
          as.data.frame.matrix(comp_tbl),
          file.path(report_tables_folder, paste0("extrastage_table_", col_name, "_", cluster_use, ".RDS")),
          COMPRESSION_FORMAT
        )
      }
    }

    tbl <- table(scrna$name, scrna@meta.data[, cluster_use])
    rowsums <- rowSums(tbl)
    tbl <- cbind(tbl, rowsums)
    colsums <- colSums(tbl)
    tbl <- rbind(tbl, colsums)
    save_object(
      tbl,
      file.path(report_tables_folder, paste0("stSample_table_cluster_",cluster_use,".RDS")),
      COMPRESSION_FORMAT
    )

    tbl <- table(scrna$stage, scrna@meta.data[, cluster_use])
    rowsums <- rowSums(tbl)
    tbl <- cbind(tbl, rowsums)
    colsums <- colSums(tbl)
    tbl <- rbind(tbl, colsums)
    save_object(
      tbl,
      file.path(report_tables_folder, paste0("stCond_table_cluster_",cluster_use,".RDS")),
      COMPRESSION_FORMAT
    )

    if(cluster_use != "singleton"){
      ## Clusters Statistics
      message("### Making cluster statistics barplot")
      df <- seutools_partition(scrna,
                               partition=fisher_cluster_name,
                               save_dir=SAVE_DIR,
                               allinone=ALLINONE)
      df <- df %>% mutate_at(vars(starts_with("pval.adjust_")),~formatC(x=.,format = "e",digits = 3))


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

      wrapp_invalid_name <- function(a_str){
              ret <- stringr::str_replace_all(a_str, " ", "____")
              ret <- stringr::str_replace_all(ret, "\t", "____")
              return(ret)
      }

      col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(df$Cluster)))
      colnames(df) <- sapply(colnames(df), wrapp_invalid_name)

      if(length(table(scrna$stage))<=2){
        nm <- scrna@tools[["meta_order"]][["stage"]][1]
        nm2 <- scrna@tools[["meta_order"]][["stage"]][2]

        aodds <- wrapp_invalid_name(glue("odds.ratio_{nm}.vs.{nm2}"))
        pval_col <- wrapp_invalid_name(glue::glue("pval.adjust_{nm}.vs.{nm2}"))

        plt <- ggplot(
          data=df, aes(x = .data[["Cluster"]], y = .data[[aodds]], fill = .data[["Cluster"]])
        ) +
        geom_bar(stat="identity") +
        coord_flip() +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_y_continuous(trans = shift_trans(1)) +
        geom_text(
          data = df,
          aes(
            x = .data[["Cluster"]],
            y = 1,
            label = .data[[pval_col]]
          ),
          position = "identity",
          size=4
        ) +
        scale_fill_manual(values=col_def) +
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

          aodds <- wrapp_invalid_name(glue("odds.ratio_{nm}.vs.{nm2}"))
          pval_col <- wrapp_invalid_name(glue("pval.adjust_{nm}.vs.{nm2}"))

          plt <- ggplot(data=df, aes(x = .data[["Cluster"]], y = .data[[aodds]], fill = .data[["Cluster"]])) +
          geom_bar(stat="identity") +
          coord_flip() +
          guides(fill = guide_legend(reverse = TRUE)) +
          scale_y_continuous(trans = shift_trans(1)) +
          geom_text(
            data = df,
            aes(
              x = .data[["Cluster"]],
              y = 1,
              label = .data[[pval_col]]
            ),
            position = "identity",
            size=4
          ) +
          scale_fill_manual(values=col_def) +
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

      ## scProportion Test
      scProportion_cluster_name <- paste0("proptest_", cluster_use)
      if(scProportion_cluster_name %in% names(scrna@tools)){
          message("### Making proptest per cluster")

          dffprop <- seutools_partition(scrna,
                                        partition=scProportion_cluster_name,
                                        save_dir=SAVE_DIR,
                                        allinone=ALLINONE)

          pairs <- combn(1:length(table(scrna$stage)), 2)
          n <- length(pairs)/2
          for (i in 1:n){
            i1 <- pairs[1:2, i][1]
            i2 <- pairs[1:2, i][2]
            nm <- scrna@tools[["meta_order"]][["stage"]][i1]
            nm2 <- scrna@tools[["meta_order"]][["stage"]][i2]
            dff <- dffprop[[glue("{nm}.vs.{nm2}")]]
            plt <- permutation_plot(dff) %+safe% ggtitle(glue("{nm} vs {nm2}, positive means more {nm2}"))
            save_ggplot_formats(
              plt=plt,
              base_plot_dir=report_plots_folder,
              plt_name=paste0("cluster_proptest_",nm,"-vs-",nm2,"_", cluster_use),
              width=9, height=7
            )
          }
      }
      ## Proportion
      message("### Making sample proportions per cluster barplot")
      name_len <- length(table(scrna@meta.data$name))
      help_sort_func <- ifelse(
        all.is.numeric(unique(scrna@meta.data[, cluster_use])),
        function(x) as.numeric(as.character(x)),
        as.character
      )
      scrna@meta.data[,cluster_use] <- help_sort_func(scrna@meta.data[,cluster_use])

      cluster_propotion_before = t(prop.table(x = table(scrna@meta.data$name, scrna@meta.data[, cluster_use]), margin = 2))
      cluster_propotion_before = cluster_propotion_before[,1:name_len]
      cluster_propotion_before_sort = cluster_propotion_before[order(
        cluster_propotion_before[,1],
        cluster_propotion_before[,2],
        decreasing=TRUE
      ),]

      col_def <- ggsci_pal(option = replicates_viridis_opt)(name_len)

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
        t(cluster_propotion_before_sort), col=col_def, xlab="cluster",
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
        t(cluster_propotion_before_sort), col=col_def, xlab="cluster",
        legend.text = colnames(cluster_propotion_before_sort),
        args.legend = list(x ='right', bty='n', inset=c(-0.13,0), xpd = TRUE),
        names.arg = rownames(cluster_propotion_before_sort),
        main = "Propotion of dataset"
      )
      dev.off()



      ## Amount Distribution
      message("### Making sample amount distribution barplot")
      #myel <- SetAllIdent(scrna, "name")

      col_def <- ggsci_pal(option = replicates_viridis_opt)(length(unique(scrna@meta.data[,"name"])))

      plt = BarPlot(
        scrna@meta.data[, cluster_use], fill = scrna$name, xlab = "Cluster",
        legend.title = "Replicate", main = "Amount of samples in each cluster")
      plt = plt + scale_fill_manual(values=col_def) + theme_minimal()
      save_ggplot_formats(
        plt=plt,
        base_plot_dir=report_plots_folder,
        plt_name=paste0("cluster_sample_distribution_",cluster_use),
        width=14, height=9
      )

    }

    ## Cell cycle phase
    message("### Making cellcycle phases umap")

    group_by <- "Phase"
    col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(scrna@meta.data[,group_by])))

    plt = DimPlot(
      scrna,
      reduction = umap_reduction,
      group.by = group_by,
      cols=col_def
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_ccycle_phase_",cluster_use),
      width=9, height=7
    )

    ## FeaturePlot
    message("### Making umap featureplot with QC elements")
    col_def <- c(base_color,pos_color)
    plt = StyleFeaturePlot(
      scrna,
      style=FEATUREPLOT_STYLE,
      features = c("percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA"),
      cols = col_def,
      order = T,
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
    plt = StyleFeaturePlot(
      scrna,
      style=FEATUREPLOT_STYLE,
      features = c("CC.Difference","G1.Score", "S.Score", "G2M.Score"),
      cols = col_def,
      order = T,
      reduction=umap_reduction,
      ncol = 2
    )
    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("umap_featureplot_ccycle_",cluster_use),
      width=9, height=7
    )

    if (DOUBLET_SWITCH != 'off'){
      message("### Making umap featureplot with doublets scores")
      plt = StyleFeaturePlot(
        scrna,
        style=FEATUREPLOT_STYLE,
        features = c("pANN"),
        cols = col_def,
        order = T,
        reduction=umap_reduction,
        ncol = 2
      )
      save_ggplot_formats(
        plt=plt,
        base_plot_dir=report_plots_folder,
        plt_name=paste0("umap_featureplot_doublets_",cluster_use),
        width=9, height=3.5
      )
    }


    ## Violin Plot
    message("### Making violinplot with QC and ccycle elements")
    group_by <- cluster_use
    if (DOUBLET_SWITCH=='off'){
       feats_to_plot <- c(
        "percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA",
        "CC.Difference", "G1.Score", "S.Score", "G2M.Score"
      )
    }else{
      feats_to_plot <- c(
        "percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA",
        "CC.Difference", "G1.Score", "S.Score", "G2M.Score", "pANN"
      )
    }
    col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(scrna@meta.data[,group_by])))
    plt = VlnPlot(
      scrna,
      features = feats_to_plot,
      group.by = group_by,
      cols = col_def,
      pt.size = 0,
      ncol = 3
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
      scrna <- safe_join_layers(scrna)  # v5: join layers before GetAssayData
      tmp_scrna <- CreateSeuratObject(counts = GetAssayDataCompat(scrna, assay = "RNA", layer = "counts"), meta.data = scrna@meta.data)
      tmp_scrna@reductions[[umap_reduction]] <- scrna@reductions[[umap_reduction]]
      group_by <- "MCA_annotate"
      col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(tmp_scrna@meta.data[,group_by])))
      plt <- DimPlot(
        object = tmp_scrna,
        reduction = umap_reduction,
        pt.size = 0.2,
        cols=col_def,
        group.by = group_by
      ) %+safe%
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
      save_object(
        plt,
        file.path(savedir, "to_plot",paste0("mca_annotate_plt_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
      save_object(
        FetchData(object = tmp_scrna, vars = c("MCA_annotate", cluster_use)),
        file.path(savedir,"to_plot", paste0("mca_annotate_info_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
    }

    ## HCL annotation
    if("HCL_annotate" %in% names(scrna@meta.data)){
      message("### Making umap with HCL annotation")
      scrna <- safe_join_layers(scrna)  # v5: join layers before GetAssayData
      tmp_scrna <- CreateSeuratObject(counts = GetAssayDataCompat(scrna, assay = "RNA", layer = "counts"), meta.data = scrna@meta.data)
      tmp_scrna@reductions[[umap_reduction]] <- scrna@reductions[[umap_reduction]]
      group_by <- "HCL_annotate"
      col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(tmp_scrna@meta.data[,group_by])))
      plt <- DimPlot(
        object = tmp_scrna,
        reduction = umap_reduction,
        pt.size = 0.2,
        cols=col_def,
        group.by = group_by
      ) %+safe%
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
      save_object(
        plt,
        file.path(savedir,"to_plot", paste0("hcl_annotate_plt_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
      save_object(
        FetchData(object = tmp_scrna, vars = c("HCL_annotate", cluster_use)),
        file.path(savedir, "to_plot", paste0("hcl_annotate_info_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
    }

    ## External Annotation
    if("external_annotation" %in% names(scrna@meta.data)){
      message("### Making umap with external annotation")
      scrna <- safe_join_layers(scrna)  # v5: join layers before GetAssayData
      tmp_scrna <- CreateSeuratObject(counts = GetAssayDataCompat(scrna, assay = "RNA", layer = "counts"), meta.data = scrna@meta.data)
      tmp_scrna@reductions[[umap_reduction]] <- scrna@reductions[[umap_reduction]]
      group_by <- "external_annotation"
      col_def <- ggsci_pal(option = cluster_viridis_opt)(length(unique(tmp_scrna@meta.data[,group_by])))
      plt <- DimPlot(
        object = tmp_scrna,
        reduction = umap_reduction,
        pt.size = 0.2,
        cols=col_def,
        group.by = group_by
      ) %+safe%
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
      save_object(
        plt,
        file.path(savedir, "to_plot",paste0("external_annotation_plt_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
      save_object(
        FetchData(object = tmp_scrna, vars = c("external_annotation", cluster_use)),
        file.path(savedir, "to_plot", paste0("external_annotation_info_",cluster_use,".RDS")),
        COMPRESSION_FORMAT
      )
    }
  }
}


