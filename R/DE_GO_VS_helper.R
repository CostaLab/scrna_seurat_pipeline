

compare_vs_loop <- function(
  cluster_use,
  all_de_list,
  all_goup_list,
  all_godown_list,
  pairs_list
){

  for(apair in pairs_list){

    tX <- apair[1]
    tY <- apair[2]
    a_vs <- glue("{tX}.vs.{tY}")

    message(sprintf("vs %s %s", tX, tY))
    de.list <- all_de_list[[a_vs]]
    de.list <- lapply(de.list, subset, subset = p_val_adj < 0.05)
    plots <- list()
    for (i in names(de.list)){

      if(nrow(de.list[[i]]) == 0){
        next
      }
      if("avg_logFC" %in% names(de.list[[i]])){ ## compatible for seurat3
        de.list[[i]]$avg_log2FC <- de.list[[i]]$avg_logFC / log(2)
      }
      x.lim <- max(abs(de.list[[i]]$avg_log2FC))
      x.lim <- c(-x.lim, x.lim)
      plots[[i]] <- GeneBarPlot(
        de.list[[i]],
        xlim = x.lim,
        main = paste0("Cluster: ", as.character(i))
      )
    }
    plots <- Filter(Negate(is.null), plots)

    if(length(plots) > 0){
      for (i in seq(1, length(plots), by = 3)){
        ni <- min(i + 2, length(plots))
        get_plts <- plots[i:ni]

        plt <- plot_grid(plotlist = get_plts, ncol = 3)
        dynamic_height <-
          2 + (0.2 * max(sapply(get_plts, function(p){nrow(p$data)})))
        # min size + (per_row_size * max_rows_in_plts_data)

        save_ggplot_formats(
          plt = plt,
          base_plot_dir = report_plots_folder,
          plt_name = glue("top10-deg-vs_{a_vs}_cluster-{cluster_use}_p{i}-{ni}"),
          width = 9,
          height = dynamic_height
        )
      }
    }

    #### {{tX}} vs {{tY}} Volcano
    de.list <- all_de_list[[a_vs]]
    help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})

    for (id in sort(help_sort_func(names(de.list)))) {
      id <- as.character(id)
      a_de <- de.list[[id]]
      if("avg_logFC" %in% names(a_de)){ ## compatible for seurat3
        a_de$avg_log2FC <- a_de$avg_logFC / log(2)
      }
      a_de$log2FC <- a_de$avg_log2FC #/ log(2)
      up <- nrow(a_de %>% filter(log2FC >= 0.25 & p_val_adj <= 0.05) )
      down <- nrow(a_de %>% filter(log2FC <= -0.25 & p_val_adj <= 0.05))
      plt <- EnhancedVolcano(
        a_de,
        x = "log2FC",
        y = "p_val_adj",
        lab = rownames(a_de),
        pCutoff = 0.05,
        FCcutoff = 0.25,
        pointSize = 1.0,
        title = glue("Volcano plot, Cluster: {id}"),
        subtitle = glue("up:{up} down:{down}")
      )

      #volcanoplot_deg_cluster-seurat_clusters_id-0.pdf
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = glue("volcanoplot_deg_vs_{a_vs}_cluster-{cluster_use}_id-{id}"),
        width = 9,
        height = 7
      )
    }

    #### {{tX}} vs {{tY}}  GO up
    go_ups <- all_goup_list[[a_vs]]
    go_ups <- go_ups[sapply(go_ups, function(x) dim(x)[1]) > 0]
    df_list <- go_ups
    df_list <- lapply(names(go_ups), function(x) go_ups[[x]]@result)
    names(df_list) <- names(go_ups)

    if(length(df_list) >= 2){

      df_list <- lapply(names(go_ups), function(x) go_ups[[x]]@result)
      for(i in 1:length(df_list)){
        assertthat::assert_that(all(rownames(df_list[[i]]) == df_list[[i]]$ID))
      }

      names(df_list) <- names(go_ups)

      pw_mtx <- pw_mtx_create(df_list)

      if(!is.null(pw_mtx)){

        col_fun <-  circlize::colorRamp2(
          c(0, 0.5, 2),
          c("purple", "black", "yellow")
        )

        plt <- Heatmap(
          pw_mtx,
          name = glue("-log10(padjust)"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = TRUE,
          col = col_fun,
          heatmap_legend_param = list(direction = "horizontal")
        )
        save_ggplot_formats(
          plt = plt,
          plt_name = glue("go_up_genes_heatmap_vs_{a_vs}_cluster-{cluster_use}"),
          units = "in",
          base_plot_dir = report_plots_folder,
          plot_obj = "ComplexHeatmap",
          width = 13,
          height = 1 + (0.3 * nrow(pw_mtx)),
          heatmap_legend_side = "top"
        )
      }
    }

    go_ups <- all_goup_list[[a_vs]]
    go_ups <- go_ups[sapply(go_ups, function(x) dim(x)[1]) > 0]

    plots <- pw_bar_plots(go_ups)

    if(length(plots) > 0){
      for (i in seq(1, length(plots), by = 2)){

        ni <- min(i + 1, length(plots))
        get_plts <- plots[i:ni]
        plt <- plot_grid(plotlist = get_plts, ncol = 2)
        dynamic_height <-
          2 + (0.3 * max(sapply(get_plts, function(p){nrow(p$data)})))

        save_ggplot_formats(
          plt = plt,
          base_plot_dir = report_plots_folder,
          plt_name = glue("top10-goup_bar_vs_{a_vs}_cluster-{cluster_use}_p{i}-{ni}"),
          width = 18,
          height = dynamic_height
        )
      }
    }

    #### {{tX}} vs {{tY}}  GO down
    go_downs <- all_godown_list[[a_vs]]
    go_downs <- go_downs[sapply(go_downs, function(x) dim(x)[1]) > 0]

    df_list <- go_downs
    df_list <- lapply(names(go_downs), function(x) go_downs[[x]]@result)
    names(df_list) <- names(go_downs)

    if(length(df_list) >= 2){

      df_list <- lapply(names(go_downs), function(x) go_downs[[x]]@result)
      for(i in 1:length(df_list)){
        assertthat::assert_that(all(rownames(df_list[[i]]) == df_list[[i]]$ID))
      }

      names(df_list) <- names(go_downs)

      pw_mtx <- pw_mtx_create(df_list)

      if(!is.null(pw_mtx)){

        col_fun <-  circlize::colorRamp2(
          c(0, 0.5, 2),
          c("purple", "black", "yellow")
        )

        plt <- Heatmap(
          pw_mtx,
          name = glue("-log10(padjust)"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = TRUE,
          col = col_fun,
          heatmap_legend_param = list(direction = "horizontal")
        )
        save_ggplot_formats(
          plt = plt,
          plt_name = glue("go_down_genes_heatmap_vs_{a_vs}_cluster-{cluster_use}"),
          units = "in",
          base_plot_dir = report_plots_folder,
          plot_obj = "ComplexHeatmap",
          width = 13,
          height = 1 + (0.3 * nrow(pw_mtx)),
          heatmap_legend_side = "top"
        )
      }
    }

    ### down genes top 10
    go_downs <- all_godown_list[[a_vs]]
    go_downs <- go_downs[sapply(go_downs, function(x) dim(x)[1]) > 0]
    plots <- pw_bar_plots(go_downs)
    if(length(plots) > 0){
      for (i in seq(1, length(plots), by = 2)){
        ni <- min(i + 1, length(plots))
        get_plts <- plots[i:ni]
        plt <- plot_grid(plotlist = get_plts, ncol = 2)
        dynamic_height <-
          2 + (0.3 * max(sapply(get_plts, function(p){nrow(p$data)})))
        save_ggplot_formats(
          plt = plt,
          base_plot_dir = report_plots_folder,
          plt_name = glue("top10-godown_bar_vs_{a_vs}_cluster-{cluster_use}_p{i}-{ni}"),
          width = 18,
          height = dynamic_height
        )

      }
    }
    #### UMAP for {{tX}} vs {{tY}}
    cluster.de.top10 <- lapply(
      de.list,
      function(x) {
        if (is.null(x[[1]])) return(NULL)
        if("avg_logFC" %in% names(x)){ ## compatible for seurat3
          x$avg_log2FC <- x$avg_logFC / log(2)
        }
        x %>% top_n(10, avg_log2FC) %>% arrange(gene)
      }
    )

    for (i in names(de.list)){
      if(nrow(cluster.de.top10[[i]]) == 0) {
        next
      }
      plt <- StyleFeaturePlot(
        scrna,
        style=FEATUREPLOT_STYLE,
        features = cluster.de.top10[[i]]$gene,
        order = TRUE,
        max.cutoff = "q95",
        cols = zero_pos_divergent_colors,
        reduction = "DEFAULT_UMAP"
      )

      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = glue("featureplot_top10_deg_vs_{a_vs}_cluster-{cluster_use}_id-{i}"),
        width = 12, height = 7
      )
    }
  }
  return(NULL)
}
