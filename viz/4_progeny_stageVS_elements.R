progeny_stageVS_elements <- function(scrna){
  cluster_use <- cluster
  progeny_cond_stage <-paste0("progeny_stage_", cluster_use)
  print(progeny_cond_stage)
  if(progeny_cond_stage %ni% names(scrna@tools)){
    stop(glue("ERORR:progeny stages comparing hasn't been calculated for cluster:{cluster_use}\n Please run [scrna_progeny_stage]!!!"))
  }

  all_progeny_list <- seutools_partition(scrna,
                                   partition=progeny_cond_stage,
                                   save_dir=SAVE_DIR,
                                   allinone=ALLINONE)

  minr <- min(sapply(all_progeny_list, function(x) min(x$r)))
  maxr <- max(sapply(all_progeny_list, function(x) max(x$r)))

  sample_stages <- scrna@tools$meta_order$stage
  list_stage = comb_list(sample_stages)


  for (apair in list_stage){
    tX <- apair[1]
    tY <- apair[2]

    a_vs <- glue("{tX}.vs.{tY}")
    progeny_df <- all_progeny_list[[a_vs]]
    help_sort_func <- ifelse(
      all.is.numeric(unique(progeny_df$CellType)),
      as.numeric,
      function(x){x}
    )
    progeny_df$CellType <- factor(
      progeny_df$CellType,
      levels= as.character(sort(unique(help_sort_func(progeny_df$CellType))))
    )

    plt <- ggplot(progeny_df, aes(y=pathway,x=CellType,fill=r)) +
            geom_tile()+
            geom_text(aes(y=pathway, x=CellType, label=tag),
                position = position_dodge(width = 0),
                hjust = 0.5, size = 7)+
            ggtitle(glue("{a_vs} r effect size")) +
            scale_fill_gradientn(colours = neg_pos_divergent_palette) +
            theme_minimal()+
            labs(x="Cluster") +
            theme(
              strip.text.x = element_text(size=28, colour="black",hjust=0),
              plot.caption = element_text(size=30, colour="black", hjust=0) #,
              # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            )

    save_ggplot_formats(
      plt=plt,
      base_plot_dir=report_plots_folder,
      plt_name=glue("progeny_heatmap_stars_vs_{a_vs}_cluster-{cluster_use}"),
      width=9, height=7
    )

    plt <- ggplot(progeny_df, aes(y=pathway,x=CellType,fill=r)) +
            geom_tile()+
            ggtitle(glue("{a_vs} r effect size")) +
            scale_fill_gradientn(colours = neg_pos_divergent_palette) +
            theme_minimal()+
            labs(x="Cluster") +
            theme(
              strip.text.x = element_text(size=28, colour="black",hjust=0),
              plot.caption = element_text(size=30, colour="black", hjust=0)#,
              # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            )
    save_ggplot_formats(plt=plt,
                    base_plot_dir=report_plots_folder,
                    plt_name=glue("progeny_heatmap_vs_{a_vs}_cluster-{cluster_use}"),
                    width=9, height=7)
  }
}

