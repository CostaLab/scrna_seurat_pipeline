external_markers_elements <- function(scrna){
  len <- length(scrna@tools$parameter)

  assertthat::assert_that(
    basename(scrna@tools$parameter[[len]]['external_file']) == basename(ext_annot_fp)
  )


  ORGAN = scrna@tools$parameter[[len]]['organ']
  SPECIES = scrna@tools$parameter[[len]]['species']

  df <- read.csv(file=ext_annot_fp, sep="\t", stringsAsFactors=F)
  df <- df[!apply(df, 1, function(x) all(x=="")), ]


  if(!(ORGAN %in% df$Tissue.of.Origin)){
       stop("exit 1 no ORGAN")
  }


  if(!(SPECIES %in% c("Human", "Mouse"))){
       stop("exit 1 NO SPECIES")
  }


  mdf = df[df$Tissue.of.Origin == ORGAN, c(glue("{SPECIES}.Gene"), "Cell.Type")]
  celltype_names <- unique(mdf$Cell.Type)
  saveRDS(
    celltype_names,
    file = file.path(report_tables_folder, "ext_annot_celltype_names.RDS")
  )

  ## External markers
  message(paste0("### ","External markers"))
  o_genes <- rownames(scrna)
  for (a_celltype in celltype_names){

    genes <- mdf[mdf$Cell.Type==a_celltype,glue::glue("{SPECIES}.Gene") ]

    message(paste0("### ",a_celltype))
    genes <- intersect(genes, o_genes)
    if (length(genes) == 0){
        next
    }
    col_def = c(base_color,pos_color)
    for (i in seq(1, length(genes), by=4)){
      ni = min(i+3, length(genes))
      p1 <- FeaturePlot(
        object = scrna,
        pt.size=0.01,
        label=T,
        label.size=2,
        features = genes[i:ni],
        reduction = "DEFAULT_UMAP",
        order = T,
        # cols = c("lightgrey", "red"),
        cols = col_def,
        ncol = 2,
        max.cutoff = 'q95'
      )
      save_ggplot_formats(
        plt=p1,
        base_plot_dir=report_plots_folder,
        plt_name=paste0("extmarkers_inte_umap_featureplot_",a_celltype,"-genes_",i,"-to-",ni),
        width=9, height=7
      )
    }


    p2 <- DotPlot(
      object= scrna,
      features = genes,
      group.by = cluster
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

    save_ggplot_formats(
      plt=p2,
      base_plot_dir=report_plots_folder,
      plt_name=paste0("extmarkers_dotplot_",a_celltype,"_groupby-",cluster),
      width=9, height=7
    )

    group_by <- cluster
    col_def <- viridis::viridis_pal(option = cluster_viridis_opt)(length(unique(scrna@meta.data[,group_by])))
    for (i in seq(1, length(genes), by=9)){
      ni = min(i+8, length(genes))
      p3 <- VlnPlot(
        object = scrna,
        pt.size=0,
        features = genes[i:ni],
        cols = col_def,
        group.by = group_by
      )
      save_ggplot_formats(
        plt=p3,
        base_plot_dir=report_plots_folder,
        plt_name=paste0("extmarkers_vlnplot_",a_celltype,"-genes_",i,"-to-",ni),
        width=9, height=7
      )
    }
  }
}

