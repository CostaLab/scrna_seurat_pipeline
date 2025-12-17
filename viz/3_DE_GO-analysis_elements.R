DE_GO_analysis_elements <- function(scrna){
  is_contigous_true_df <- function(is_sigs){
    ret_df <- data.frame(keep = FALSE, avgIdx = -1)
    if(any(is_sigs) & table(is_sigs)['TRUE'] == 1){
      ret_df$keep = TRUE
      ret_df$avgIdx = which(is_sigs == TRUE)
      return(ret_df)
    }
    return(ret_df)
  }

  de_cluster_name <- paste0("de_", cluster)
  go_cluster_name <- paste0("go_", cluster)
  help_sort_func <- ifelse(
    all.is.numeric(unique(scrna@meta.data[, cluster])),
    function(x) as.numeric(as.character(x)),
    as.character
  )
  scrna@meta.data[, cluster] <- help_sort_func(scrna@meta.data[, cluster])
  col_def <- rev(ggsci_pal(option = cluster_viridis_opt)(length(unique(scrna@meta.data[, cluster]))))


  hallmark_cluster_name <- ""
  kegg_cluster_name <- ""
  reactome_cluster_name <- ""

  if("DEGO" %in% funcs){
    if(de_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:DE hasn't been calculated for cluster:{cluster}\n Please run [scrna_markergenes]!!!"))
    }
    if(go_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:GO hasn't been calculated for cluster:{cluster}\n Please run [scrna_go]!!!"))
    }
  }
  if("hallmark" %in% funcs){
    hallmark_cluster_name <- paste0("hallmark_", cluster)
    if(hallmark_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:hallmark hasn't been calculated for cluster {cluster}\n Please run [scrna_hallmark]!!!"))
    }
  }
  if("KEGG" %in% funcs){
    kegg_cluster_name <- paste0("kegg_", cluster)
    if(kegg_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:kegg hasn't been calculated for cluster {cluster}\n Please run [scrna_kegg]!!!"))
    }
  }
  if("Reactome" %in% funcs){
    reactome_cluster_name <- paste0("reactome_", cluster)
    if(reactome_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:reactome hasn't been calculated for cluster {cluster}\n Please run [scrna_reactome]!!!"))
    }
  }
  if("progeny" %in%funcs){
    progeny_cluster_name <- paste0("progeny_", cluster)
    if(progeny_cluster_name %ni% names(scrna@tools)){
      stop(glue("ERROR:progeny hasn't been calculated for cluster {cluster}\n Please run [scrna_progeny]!!!"))
    }
    progeny_df <- seutools_partition(scrna,
                                     partition=progeny_cluster_name,
                                     save_dir=SAVE_DIR,
                                     allinone=ALLINONE)

    help_sort_func <- ifelse(
      all.is.numeric(unique(progeny_df$CellType)),
      function(x) as.numeric(as.character(x)),
      function(x){x}
    )

    progeny_df$CellType <- factor(
      progeny_df$CellType,
      levels = as.character(
        sort(unique(help_sort_func(progeny_df$CellType)))
      )
    )

    plt <-
      ggplot(progeny_df, aes(y = pathway, x = CellType, fill = r)) +
      geom_tile() +
      scale_fill_gradientn(colours = neg_pos_divergent_palette) +
      ggtitle(glue("{cluster} r effect size")) +
      labs(x = "Cluster") +
      # scale_fill_distiller(palette ="RdBu", direction = -1) +
      theme_minimal() +
      theme(
        strip.text.x = element_text(size = 28, colour = "black", hjust = 0),
        plot.caption = element_text(size = 30, colour = "black", hjust = 0)#,
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      )

    save_ggplot_formats(
      plt = plt,
      base_plot_dir = report_plots_folder,
      plt_name = paste0("progeny_r_effect_heatmap-", cluster),
      width = 9, height = 7
    )


    plt <-
      ggplot(progeny_df, aes(y = pathway, x = CellType, fill = r)) +
      geom_tile() +
      geom_text(
        aes(y = pathway, x = CellType, label = tag),
        position = position_dodge(width = 0),
        hjust = 0.5, size = 4
      ) +
      scale_fill_gradientn(colours = neg_pos_divergent_palette) +
      ggtitle(glue("{cluster} r effect size")) +
      labs(x = "Cluster") +
      # scale_fill_distiller(palette ="RdBu", direction = -1) +
      theme_minimal() +
      theme(
        strip.text.x = element_text(size = 28, colour = "black", hjust = 0),
        plot.caption = element_text(size = 30, colour = "black", hjust = 0)#,
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      )

    save_ggplot_formats(
      plt = plt,
      base_plot_dir = report_plots_folder,
      plt_name = paste0("progeny_r_effect_text_heatmap-", cluster),
      width = 9, height = 7
    )
  }



  if("DEGO" %in% funcs){
    # DE plots
    cluster_de <- seutools_partition(scrna,
                                     partition=de_cluster_name,
                                     save_dir=SAVE_DIR,
                                     allinone=ALLINONE)


    cluster_de <- cluster_de[sapply(cluster_de, function(m) nrow(m) > 0)]

    cluster_de_top10 <- lapply(cluster_de, function(x) {
        if("avg_logFC" %in% names(x)){ ## compatible for seurat3
          x$avg_log2FC <- x$avg_logFC / log(2)
        }
        x %>% top_n(10, avg_log2FC) %>% arrange(-avg_log2FC)
    })


    ## top10 DE heatmaps
    genes <- as.vector(unlist(sapply(cluster_de_top10, function(x)x$gene)))
    scrna <- ScaleData(scrna, rownames(scrna))


    scrna_to_plot <- function(scrna){
      if(ncol(scrna) <= 20000){
        return(scrna)
      }
      return(subset(scrna, cells = sample(colnames(scrna), size = 20000)))
    }

    plthm = DoHeatmap(
      ## >30k will failed to plot, here we subset when cell number > 20,000
      object = scrna_to_plot(scrna),
      features = genes,
      group.by = cluster,
      group.colors = col_def,
      disp.min = -2,
      disp.max = 2,
      slot = "scale.data",
      assay = "RNA",
      raster = FALSE,
      combine = TRUE
    ) +
    ggtitle("Marker genes for each cluster") +
    NoLegend()

    save_ggplot_formats(
      plt = plthm,
      base_plot_dir = report_plots_folder,
      plt_name = paste0("heatmap_top10-de-genes_cluster-", cluster),
      width = 13,
      height = 12
    )


    # GeneBarPlots
    ## Plot the top 10 DE genes in each cluster.
    plots = list()

    help_sort_func <- ifelse(
      all.is.numeric(names(cluster_de)),
      as.numeric,
      function(x){x}
    )

    for (id in sort(help_sort_func(names(cluster_de)))) {
      id <- as.character(id)
      cluster_genes <- cluster_de_top10[[id]]
      if("avg_logFC" %in% names(cluster_de[[id]])){ ## compatible for seurat3
        cluster_de[[id]]$avg_log2FC <- cluster_de[[id]]$avg_logFC / log(2)
      }
      x_lim = max(abs(cluster_de[[id]]$avg_log2FC))
      x_lim <- c(-x_lim, x_lim)
      plots[[id]] <- GeneBarPlot(
        cluster_de[[id]],
        xlim = x_lim,
        main = paste0("Cluster: ", id)
      )
    }

    if(length(plots) > 0){
      for (i in seq(1, length(plots), by = 4)){
        ni <- min(i + 3, length(plots))
        get_plts <- plots[i:ni]
        plt <- plot_grid(plotlist = get_plts, ncol = 4)

        dynamic_height <-
          2 + (0.2 * max(sapply(get_plts, function(p){nrow(p$data)})))

        save_ggplot_formats(
          plt = plt,
          base_plot_dir = report_plots_folder,
          plt_name = paste0(
            "top10-de-genes_per_cluster-", cluster, "_p", i, "-", ni
          ),
          width = 13,
          height = dynamic_height
        )
      }
    }

    # volcano plots
    help_sort_func <- ifelse(
      all.is.numeric(names(cluster_de)),
      as.numeric,
      function(x){x}
    )

    for (id in sort(help_sort_func(names(cluster_de)))) {
      id = as.character(id)
      a_de <- cluster_de[[id]]
      if("avg_logFC" %in% names(a_de)){ ## compatible for seurat3
        a_de$avg_log2FC <- a_de$avg_logFC / log(2)
      }
      a_de$log2FC <- a_de$avg_log2FC # / log(2)
      up <- nrow(a_de %>% filter(log2FC >= 1 & p_val_adj <= 0.05) )
      down <- nrow(a_de %>% filter(log2FC <= -1 & p_val_adj <= 0.05))
      plt <- EnhancedVolcano(
        a_de,
        x = "log2FC",
        y = "p_val_adj",
        lab = a_de$gene,
        pointSize = 1.0,
        pCutoff = 0.05,
        title = glue("Volcano plot, Cluster: {id}"),
        subtitle = glue("up:{up} down:{down}")
      )
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = paste0(
          "volcanoplot_deg_cluster-",
          cluster, "_id-",
          URLencode_escape(id)
        ),
        width = 9,
        height = 7
      )
    }

    ## DE genes on UMAP plot
    Idents(scrna) <- cluster
    for (i in names(cluster_de)){
      #plots = list()
      print(sprintf("Cluster: %s", i))
      ps <- StyleFeaturePlot(
        scrna,
        features = cluster_de_top10[[as.character(i)]]$gene,
        label = TRUE,
        max.cutoff = "q95",
        order = TRUE,
        label.size = 2,
        cols = zero_pos_divergent_colors,
        reduction = "DEFAULT_UMAP"
      )
      save_ggplot_formats(
        plt = ps,
        base_plot_dir = report_plots_folder,
        plt_name = paste0(
          "umap_featureplot_top10deg_cluster-",
          cluster, "_id-",
          URLencode_escape(i)
        ),
        width = 12, height = 7
      )
    }
  }

  #### Genesets
  # FIXME if genesets not requested, this block shouldn't run
  if("Genesets" %in% funcs){
    for(nm in scrna@tools$genesets){
      message("nm: ", nm, "    ", appendLF = FALSE)
      if(nm %ni% names(scrna@meta.data)){
        next
      }
      df <- data.frame(
        nm = scrna@meta.data[, nm],
        Cluster = as.character(scrna@meta.data[, cluster]),
        stringsAsFactors = FALSE)


      df.s <- melt(df, id.vars = c("Cluster"))
      df.s[df.s == -Inf] <- 0

      min_x <- min(df.s$value)
      max_x <- max(df.s$value)

      ## ridges
      plt <- ggplot(
          df.s,
          aes(
            x = value,
            y = Cluster,
            color = Cluster,
            point_color = Cluster,
            fill = Cluster
          )
        ) +
        geom_density_ridges(
          jittered_points = FALSE,
          scale = .95,
          rel_min_height = .01,
          alpha = 0.5
        ) +
        scale_y_discrete(expand = c(.01, 0)) +
        scale_fill_manual(values = col_def) +
        scale_color_manual(values = col_def, guide = "none") +
        scale_discrete_manual("point_color", values = col_def, guide = "none") +
        theme_ridges(center = TRUE) +
        xlim(min_x, max_x) + ylab("") + ggtitle(glue("{nm}")) +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = paste0("Genesets_ridges_", nm, "-", cluster),
        width = 9, height = 7
      )


      ## Vln
      plt <- VlnPlot(
        scrna, features = nm,
        pt.size = 0,
        group.by = cluster,
        cols = col_def) + ggtitle(nm)

      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = paste0("Genesets_violin_", nm, "-", cluster),
        width = 9, height = 7
      )


    ## feature
      plt <- StyleFeaturePlot(
        scrna, features = nm,
        reduction = "DEFAULT_UMAP",
        max.cutoff = "q95",
        order = TRUE,
        cols = zero_pos_divergent_colors
      ) + ggtitle(nm)
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = report_plots_folder,
        plt_name = paste0("Genesets_feature_", nm, "-", cluster),
        width = 9, height = 7
      )

    }
  }

  #############################################################
  ## Term enrichment analysis (GO, hallmark, KEGG, Reactome) ##
  #############################################################
  enrch_analysis_vec = c()

  if("DEGO" %in% funcs)
    enrch_analysis_vec = c(enrch_analysis_vec, go_cluster_name)
  if("hallmark" %in% funcs)
    enrch_analysis_vec = c(enrch_analysis_vec, hallmark_cluster_name)
  if("KEGG" %in% funcs)
    enrch_analysis_vec = c(enrch_analysis_vec, kegg_cluster_name)
  if("Reactome" %in% funcs)
    enrch_analysis_vec = c(enrch_analysis_vec, reactome_cluster_name)

  for(enrich_cluster_name in enrch_analysis_vec){

    enrich_type <- gsub(paste0("_", cluster, "$"), "", enrich_cluster_name)
    ## TERM up analysis
    term_up_list <-
      if(enrich_cluster_name == go_cluster_name){
        seutools_partition(scrna,
                           partition=go_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$goup
      }else if(enrich_cluster_name == hallmark_cluster_name){
        seutools_partition(scrna,
                           partition=hallmark_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$hallmarkup
      }else if(enrich_cluster_name == kegg_cluster_name){
        seutools_partition(scrna,
                           partition=kegg_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$keggup
      }else if(enrich_cluster_name == reactome_cluster_name){
        seutools_partition(scrna,
                           partition=reactome_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$reactomeup
      }


    ## TERM down analysis
    term_down_list <-
      if(enrich_cluster_name == go_cluster_name){
        seutools_partition(scrna,
                           partition=go_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$godown
      }else if(enrich_cluster_name == hallmark_cluster_name){
        seutools_partition(scrna,
                           partition=hallmark_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$hallmarkdown
      }else if(enrich_cluster_name == kegg_cluster_name){
        seutools_partition(scrna,
                           partition=kegg_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$keggdown
      }else if(enrich_cluster_name == reactome_cluster_name){
        seutools_partition(scrna,
                           partition=reactome_cluster_name,
                           save_dir=SAVE_DIR,
                           allinone=ALLINONE)$reactomedown
      }

    for(term_direction in c("up", "down")){

      term_direction_list = switch(
        term_direction,
        "up" = term_up_list,
        "down" = term_down_list
      )

      if(length(term_direction_list)  == 0){next}

      df_list <- lapply(
        1:length(term_direction_list),
        function(x) term_direction_list[[x]]@result
      )
      names(df_list) <- names(term_direction_list)

      for(i in 1:length(df_list)){
          assertthat::assert_that(
            all(rownames(df_list[[i]]) == df_list[[i]]$ID)
          )
      }
      union_TermID <- Reduce(union, lapply(df_list, function(x) x$ID))

      union_df <- do.call(rbind, df_list)
      union_df$GeneRatio <- 0
      union_df$BgRatio <- 0
      union_df$pvalue <- 1
      union_df$p.adjust <- 1
      union_df$qvalue <- 1
      union_df$geneID <- ""
      union_df$Count <- 0

      union_df <- union_df[!duplicated(union_df$ID), ]
      rownames(union_df) <- union_df$ID
      for(x in 1:length(df_list)){
          rest_ids <- setdiff(rownames(union_df), rownames(df_list[[x]]))
          df_list[[x]] <- rbind(df_list[[x]], union_df[rest_ids, ])
      }
      #rest_ids[1:10]
      filtered_term <- c()

      avgIdx <- list()
      for(TermID in  union_TermID){
          is_sigs <- sapply(df_list, function(x)x[x$ID == TermID,]$p.adjust < 0.05)
          is_true_df <- is_contigous_true_df(is_sigs)
          if(is_true_df$keep){
              filtered_term <- c(filtered_term, TermID)
              avgIdx[[ union_df[union_df$ID == TermID,]$Description ]] <- is_true_df$avgIdx
          }

      }
      if(length(filtered_term) > 5){
        df_list <- lapply(df_list, function(x) x %>% filter(ID %in% filtered_term) )
        nms <- names(df_list)
        df_list <- lapply(names(df_list), function(x) df_list[[x]] %>% mutate(name = x))
        names(df_list) <- nms
        df_list_select <- lapply(
          1:length(df_list),
          function(x){
            df_list[[x]] %>%
              filter(p.adjust < 0.05) %>%
              top_n(wt = -log10(p.adjust), n = 5) %>%
              arrange(log10(p.adjust))
          }
        )
        df_list_select <- lapply(
          df_list_select, function(x)x[1:min(5, nrow(x)), ]
        )

        all_names <- as.vector(
          unlist(
            sapply(
              1:length(df_list_select),
              function(x){df_list_select[[x]]$ID}
            )
          )
        )

        pdf_list <- lapply(
          1:length(df_list), function(x) subset(df_list[[x]], ID %in%all_names)
        )

        mdf <- do.call(rbind, pdf_list)
        pmdf <- mdf[, c("Description", "name", "p.adjust")]
        pmdf$name <- factor(pmdf$name, levels = names(df_list))

        pmtx <- reshape2::dcast(pmdf,  Description ~ name, value.var = "p.adjust")

        rownames(pmtx) <- pmtx$Description
        pmtx$Description <- NULL
        help_mtx <- pmtx
        help_mtx[help_mtx >= 0.05] = 1000
        help_mtx[help_mtx < 0.05] = 1
        help_mtx <- help_mtx[do.call(order, help_mtx),]
        #matrix_list[[pw]] <- pmtx[rownames(help_mtx), ]
        pmtx <- -log10(pmtx)
        pmtx[pmtx>2] = 2
        pmtx <- pmtx[rownames(help_mtx), ]
        col_fun <-  circlize::colorRamp2(
          c(0, 0.5, 2),
          c("purple", "black", "yellow")
        )
        plthm <- Heatmap(
          as.matrix(pmtx)[order(unlist(avgIdx[rownames(help_mtx)])), ],
          name = glue("-log10(padjust)"),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = TRUE,
          col = col_fun,
          heatmap_legend_param = list(direction = "horizontal")
        )

        save_ggplot_formats(
          plt = plthm,
          units = "in",
          base_plot_dir = report_plots_folder,
          plt_name = paste0(
            "heatmap_",enrich_type,"-",
            term_direction,"_cluster-", cluster
          ),
          plot_obj = "ComplexHeatmap",
          width = 10,
          height = 1 + (0.3 * nrow(pmtx)),
          heatmap_legend_side  =  "top"
        )
      }

      ### direction genes top 10
      term_plot_list = lapply(
        term_direction_list,
        function(x){
          df = x@result
          if (is.list(df) && length(df) == 0){
            log_m = as.data.frame(list())
            return(log_m)
          }
          log_m = as.data.frame(-log10(df$p.adjust))
          log_m$names = as.factor(sapply(df$Description, function(y){
              y <- as.character(trimws(y))
              return(y) }))
          log_m$show_names = as.factor(sapply(df$Description, function(y){
              y <- as.character(trimws(y))
              y <- ifelse(nchar(y) <= 33,  y, paste0(substr(y, 1, 30), "..."))
              return(y) }))
          log_m <- log_m[order(log_m[,1], decreasing = TRUE),]
          showCatetermry = min(length(log_m[,1]), 10)
          log_m <- log_m[1:showCatetermry, ]
          log_m <- log_m[order(log_m[,1], decreasing = FALSE),]
          return(log_m)
        }
      )

      ### direction genes plot
      plots <- lapply(
        seq_along(term_plot_list),
        function(y, i) {
          term_df <- y[[i]]
          if(length(col) == 0)
            return(NULL)
          ggplot(
            term_df,
            aes(reorder(x = term_df[,2], term_df[,1]), y = term_df[,1])
          ) +
          geom_bar(
            stat = "identity",
            fill = ifelse(term_direction == "up", pos_color, neg_color)
          ) +
          geom_hline(
            yintercept = -log10(0.05),
            linetype = 2, size = 1, color = "grey50"
          ) +
          ggtitle(
            paste0(enrich_type, " ", term_direction, ", Cluster:", names(y)[i])
          ) +
          theme_minimal(base_size = 16) +
          theme(
            axis.text.y  = element_text(size = 20),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank()
          ) +
          labs(y = "") +
          scale_y_continuous(name = "-log10(p-value)") +
          scale_x_discrete(breaks = term_df[,2], labels = term_df[,3]) +
          coord_flip()
        },
        y = term_plot_list
      )

      #plt = plot_grid(plotlist = plots, ncol=2)
      if(length(plots) > 0){
          for (i in seq(1, length(plots), by = 2)){
            ni <- min(i + 1, length(plots))
            get_plts <- plots[i:ni]

            plt <- plot_grid(plotlist = get_plts, ncol = 2)

            dynamic_height <-
              2 + (0.2 * max(sapply(get_plts, function(p){nrow(p$data)})))

            save_ggplot_formats(
              plt = plt,
              base_plot_dir = report_plots_folder,
              plt_name = paste0(
                enrich_type, "_",
                term_direction, "_genes_barplot_cluster-",
                cluster, "_p", i, "-", ni
              ),
              width = 16,
              height = dynamic_height
            )
          }
      }
    }
  }
}

