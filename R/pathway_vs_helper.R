pathway_vs_loop <- function(
  cluster_use,
  pw_param,
  pathways,
  pairs_list
){

  for(pw in pathways){
    message("pathway: ", pw)
    for(direction in c("up", "down")){

      pw_partition = seutools_partition(scrna, pw_param, SAVEDIR, allinone=scnra@tools$allinone)
      all_pw_direction_list <- pw_partition[[glue("{pw}{direction}")]]

      for(apair in pairs_list){
        tX <- apair[1]
        tY <- apair[2]
        a_vs <- glue("{tX}.vs.{tY}")

        pw_directions <- all_pw_direction_list[[a_vs]]
        pw_directions <- pw_directions[sapply(pw_directions, function(x) dim(x)[1]) > 0]

        df_list <- pw_directions
        df_list <- lapply(names(pw_directions), function(x) pw_directions[[x]]@result)
        names(df_list) <- names(pw_directions)

        if(length(df_list) >= 2){

          df_list <- lapply(names(pw_directions), function(x) pw_directions[[x]]@result)
          for(i in 1:length(df_list)){
            assertthat::assert_that(all(rownames(df_list[[i]]) == df_list[[i]]$ID))
          }

          names(df_list) <- names(pw_directions)
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
              units = "in",
              base_plot_dir = report_plots_folder,
              plt_name=glue("{pw}_{direction}_genes_heatmap_vs_{a_vs}_cluster-{cluster_use}"),
              plot_obj = "ComplexHeatmap",
              width = 13,
              height = 1 + (0.3 * nrow(pw_mtx)),
              heatmap_legend_side = "top"
            )
          }
        }

        pw_directions <- all_pw_direction_list[[a_vs]]
        pw_directions <- pw_directions[sapply(pw_directions, function(x) dim(x)[1]) > 0]

        plots <- pw_bar_plots(pw_directions)

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
              plt_name = glue("top10-{pw}{direction}_bar_vs_{a_vs}_cluster-{cluster_use}_p{i}-{ni}"),
              width = 16,
              height = dynamic_height
            )

          }
        }
      }
    }
  }
  return(NULL)
}

