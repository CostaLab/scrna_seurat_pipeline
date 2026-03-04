DEGO_stageVS_elements <- function(scrna){

  cluster_use <- cluster
  dego_stage_name <- paste0("dego_stage_", cluster_use)

  if(dego_stage_name %ni% names(scrna@tools)){
    stop(glue("ERROR: DE&GO stages comparing hasn't been calculated for cluster:{cluster_use}\n Please run [scrna_dego_name]!!!"))
  }

  dego_stage_df <- seutools_partition(scrna,
                                      partition=dego_stage_name,
                                      save_dir=SAVE_DIR,
                                      allinone=ALLINONE)

  all_de_list          <- dego_stage_df$de
  all_goup_list        <- dego_stage_df$goup
  all_godown_list      <- dego_stage_df$godown

  stage_names <- scrna@tools$meta_order$stage
  list_stages <- comb_list(stage_names)

  compare_vs_loop(
    cluster_use = cluster_use,
    all_de_list = all_de_list,
    all_goup_list = all_goup_list,
    all_godown_list = all_godown_list,
    pairs_list = list_stages
  )
}

DEGO_extrastageVS_elements <- function(scrna) {
  cluster_use <- cluster
  esc <- scrna@tools[["extra_stage_cols"]]
  if (is.null(esc) || length(esc) == 0) return(NULL)

  for (col_name in names(esc)) {
    tool_key <- sprintf("dego_extrastage_%s_%s", col_name, cluster_use)
    if (tool_key %ni% names(scrna@tools)) {
      message(sprintf("Skipping extra stage viz for %s: not computed", col_name))
      next
    }
    dego_df <- seutools_partition(scrna,
                                  partition = tool_key,
                                  save_dir = SAVE_DIR,
                                  allinone = ALLINONE)
    defn <- esc[[col_name]]
    meta_col <- if (is.list(defn)) col_name else defn
    groups <- na.omit(unique(as.character(scrna@meta.data[, meta_col])))
    pairs_list <- comb_list(groups)

    compare_vs_loop(
      cluster_use = cluster_use,
      all_de_list = dego_df$de,
      all_goup_list = dego_df$goup,
      all_godown_list = dego_df$godown,
      pairs_list = pairs_list,
      prefix = paste0("extrastage_", col_name, "_")
    )
  }
}

