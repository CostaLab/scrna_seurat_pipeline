pathway_extrastageVS_elements <- function(scrna){
  cluster_use <- cluster

  dic <- c(
    "hallmark_extrastage" = "hallmark",
    "reactome_extrastage" = "reactome",
    "kegg_extrastage" = "kegg"
  )

  pathways <- dic[EXEC_PLAN]
  pathways <- pathways[!is.na(pathways)]

  esc <- scrna@tools[["extra_stage_cols"]]
  if (is.null(esc) || length(esc) == 0) return(NULL)

  for (col_name in names(esc)) {
    pathway_tool_key <- sprintf("pathway_extrastage_%s_%s", col_name, cluster_use)
    if (pathway_tool_key %ni% names(scrna@tools)) {
      message(sprintf("Skipping pathway extra stage viz for %s: not computed", col_name))
      next
    }

    defn <- esc[[col_name]]
    meta_col <- if (is.list(defn)) col_name else defn
    groups <- sort(na.omit(unique(as.character(scrna@meta.data[, meta_col]))))
    pairs_list <- comb_list(groups)

    pathway_vs_loop(
      scrna = scrna,
      cluster_use = cluster_use,
      pw_param = pathway_tool_key,
      pathways = pathways,
      pairs_list = pairs_list,
      save_dir = SAVE_DIR,
      prefix = glue("extrastage_{col_name}_")
    )
  }
}