pathway_stage_elements <- function(scrna){
  cluster_use <- cluster

  dic <- c("hallmark_stage"="hallmark",
           "reactome_stage"="reactome",
           "kegg_stage"="kegg")

  pathways <- dic[EXEC_PLAN]
  pathways <-pathways[!is.na(pathways)]


  for(pw in pathways ){
    message("pathway: ", pw)
    for(direction in c("up", "down")){
      pw_sample_stage <-paste0("pathway_stage_", cluster_use)
      if(pw_sample_stage %ni% names(scrna@tools)){
        stop(glue("ERORR:{pw} samples comparing hasn't been calculated for cluster:{cluster_use}\n
                  Please run [scrna_pathway_stage]!!!"))
      }
      all_pw_direction_list <-  scrna@tools[[pw_sample_stage]][[glue("{pw}{direction}")]]

      sample_stages <- scrna@tools$meta_order$stage
      list_stage = comb_list(sample_stages)
      for(apair in list_stage){
        tX = apair[1]
        tY = apair[2]
        a_vs <- glue("{tX}.vs.{tY}")

        pw_directions <- all_pw_direction_list[[a_vs]]
        pw_directions <- pw_directions[sapply(pw_directions, function(x) dim(x)[1]) > 0]

        df_list <- pw_directions
        df_list <- lapply(names(pw_directions), function(x) pw_directions[[x]]@result)
        names(df_list) <- names(pw_directions)

        if(length(df_list) >= 2){


          df_list <- lapply(names(pw_directions), function(x) pw_directions[[x]]@result)
          for(i in 1:length(df_list)){
            assertthat::assert_that( all(rownames(df_list[[i]]) == df_list[[i]]$ID))
          }


          names(df_list) <- names(pw_directions)
          pw_mtx <- pw_mtx_create(df_list)
          if(!(is.null(pw_mtx))){
            col_fun <-  circlize::colorRamp2(c(0, 0.5, 2), c("purple", "black", "yellow"))
            plt <- Heatmap(pw_mtx,
                           name = glue("-log10(padjust)"),
                           cluster_columns = F,
                           cluster_rows = F,
                           show_row_names=T,
                           col=col_fun,
                           heatmap_legend_param = list(direction = "horizontal")
            )
            save_ggplot_formats(plt=plt,
                                units="px",
                                base_plot_dir=report_plots_folder,
                                plt_name=glue("{pw}_{direction}_genes_heatmap_vs_{a_vs}_cluster-{cluster_use}"),
                                plot_obj="ComplexHeatmap",
                                width=1200, height=1800, heatmap_legend_side = "top")
          }

        }

        pw_directions <- all_pw_direction_list[[a_vs]]
        pw_directions <- pw_directions[sapply(pw_directions, function(x) dim(x)[1]) > 0]

        plots <- pw_bar_plots(pw_directions)
        if(length(plots) > 0){
          for (i in seq(1, length(plots), by=2)){
            ni = min(i+1, length(plots))
            plt <-plot_grid(plotlist=plots[i:ni], ncol=2)
            save_ggplot_formats(plt=plt,
                                base_plot_dir=report_plots_folder,
                                plt_name=glue("top10-{pw}{direction}_bar_vs_{a_vs}_cluster-{cluster_use}_p{i}-{ni}"),
                                width=18, height=10)

          }
        }
      }
    }
  }
}
