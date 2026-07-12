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
  save_object(
    object = celltype_names,
    file_name = file.path(report_tables_folder, "ext_annot_celltype_names.RDS"),
    file_format = COMPRESSION_FORMAT
  )

  ## External markers (Seurat v5: MAGIC_RNA is optional; use RNA if not present)
  message(paste0("### ","External markers"))
  scrna <- safe_join_layers(scrna)
  Idents(scrna) <- cluster
  assay_use <- "RNA"
  if ("MAGIC_RNA" %in% names(scrna@assays)) {
    assay_use <- "MAGIC_RNA"
  } else if (!ALLINONE && "MAGIC_RNA" %in% names(scrna@tools$assay_info)) {
    scrna <- seu_assay(scrna, assay = "MAGIC_RNA", SAVE_DIR, allinone = FALSE)
    assay_use <- "MAGIC_RNA"
  }
  DefaultAssay(scrna) <- assay_use
  o_genes <- rownames(scrna)
  for (a_celltype in celltype_names){

    genes <- mdf[mdf$Cell.Type==a_celltype,glue::glue("{SPECIES}.Gene") ]

    message(paste0("### ",a_celltype))
    genes <- intersect(genes, o_genes)
    genes_filter = c()
    for(gene in genes){
      if(sum(GetAssayDataCompat(scrna, assay = "RNA", layer = "counts")[gene, ]) > 0){
        genes_filter = c(genes_filter, TRUE)
      }else{
        genes_filter = c(genes_filter, FALSE)
      }
    }
    genes = genes[genes_filter]

    if (length(genes) == 0){
        next
    }
    
    for (i in seq(1, length(genes), by = 4)){
      ni = min(i + 3, length(genes))
      p1 <- StyleFeaturePlot(
        object = scrna,
        style = FEATUREPLOT_STYLE,
        pt.size = 0.1,
        label = TRUE,
        label.size = 2,
        features = genes[i:ni],
        reduction = "DEFAULT_UMAP",
        order = TRUE,
        cols = zero_pos_divergent_colors,
        max.cutoff = "q95",
        ncol = 2
      )
      save_ggplot_formats(
        plt = p1,
        base_plot_dir = report_plots_folder,
        plt_name = paste0(
          "extmarkers_inte_umap_featureplot_",
          a_celltype, "-genes_", i, "-to-", ni
        ),
        width=9, height=7
      )
    }


    p2 <- DotPlot(
      object = scrna,
      features = genes,
      group.by = cluster,
      cols = zero_pos_divergent_colors
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

    save_ggplot_formats(
      plt = p2,
      base_plot_dir = report_plots_folder,
      plt_name = paste0(
        "extmarkers_dotplot_", a_celltype, "_groupby-", cluster
      ),
      width = 9, height = 7
    )

    group_by <- cluster
    col_def <- rev(ggsci_pal(option = cluster_viridis_opt)(
        length(unique(scrna@meta.data[,group_by]))
    ))
    for (i in seq(1, length(genes), by = 9)){
      ni = min(i + 8, length(genes))
      p3 <- VlnPlot(
        object = scrna,
        pt.size = 0,
        features = genes[i:ni],
        cols = col_def,
        group.by = group_by
      )
      save_ggplot_formats(
        plt = p3,
        base_plot_dir = report_plots_folder,
        plt_name = paste0(
          "extmarkers_vlnplot_", a_celltype, "-genes_", i, "-to-", ni
        ),
        width = 9, height = 7
      )
    }
  }
}

