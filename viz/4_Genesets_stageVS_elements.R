Genesets_stageVS_elements <- function(scrna){
  cluster_use <- cluster
  ## Ridges plots
  for(nm in scrna@tools$genesets){
    if(nm %ni% names(scrna@meta.data)){
      next
    }
    message("nm: ", nm, "   ", appendLF=F)
    df <- data.frame(
                     nm = scrna@meta.data[, nm],
                     Cluster = as.character(scrna@meta.data[, cluster_use]),
                     stage = as.character(scrna@meta.data[, "stage"]),
                     stringsAsFactors=F)
    df.s <- melt(df, id.vars = c("Cluster", "stage"))
    df.s[df.s == -Inf] <- 0

    min_x <- min(df.s$value)
    max_x <- max(df.s$value)
    col_def <- ggsci_pal(option = "futurama")(length(unique(scrna@meta.data[, "stage"])))
    plt <- ggplot(df.s, aes(x=value, y=Cluster, color=stage, point_color=stage, fill=stage)) +
      geom_density_ridges(jittered_points=FALSE, scale = .95, rel_min_height = .01, alpha=0.5) +
      scale_y_discrete(expand = c(.01, 0)) +
      scale_fill_manual(values = col_def) +
      scale_color_manual(values = col_def, guide = "none") +
      scale_discrete_manual("point_color", values = col_def, guide = "none") +
      guides(fill = guide_legend(override.aes = list(fill = col_def,
             color = NA, point_color = NA))) +
      theme_ridges(center = TRUE) +
      xlim(min_x, max_x) + ylab("") + ggtitle(glue("{nm}")) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))

      save_ggplot_formats(plt=plt,
                    base_plot_dir=report_plots_folder,
                    plt_name=glue("Geneset_ridges_vs_1v1_{nm}_cluster-{cluster_use}"),
                    width=9, height=7)
  }
  sample_stages <- scrna@tools$meta_order$stage
  list_1v1 = comb_list(sample_stages)
  scrna$stage <- as.character(scrna$stage)
  for(apair in list_1v1){
    tX = apair[1]
    tY = apair[2]
    a_vs <- glue("{tX}.vs.{tY}")
    Idents(scrna) <- "stage"
    for(nm in scrna@tools$genesets){
      plt <- VlnPlot(scrna, idents=c(tX, tY),
                     pt.size=0, features=nm,
                     cols = col_def,
                     split.plot = TRUE,
                     split.by = "stage",
                     group.by=cluster_use)
      save_ggplot_formats(plt=plt,
                          base_plot_dir=report_plots_folder,
                          plt_name=glue("Geneset_violin_vs_{a_vs}_{nm}_cluster-{cluster_use}"),
                          width=9, height=7)
    }
  }
}
