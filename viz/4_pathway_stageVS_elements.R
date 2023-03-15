pathway_stage_elements <- function(scrna){
  cluster_use <- cluster

  dic <- c(
    "hallmark_stage" = "hallmark",
    "reactome_stage" = "reactome",
    "kegg_stage" = "kegg"
  )

  pathways <- dic[EXEC_PLAN]
  pathways <-pathways[!is.na(pathways)]

  pw_sample_stage <- paste0("pathway_stage_", cluster_use)

  if(pw_sample_stage %ni% names(scrna@tools)){
    stop(glue("ERROR:{pw} samples comparing hasn't been calculated for cluster:{cluster_use}\n
              Please run [scrna_pathway_stage]!!!"))
  }

  sample_stages <- scrna@tools$meta_order$stage
  list_stage <- comb_list(sample_stages)

  pathway_vs_loop(
    scrna = scrna,
    cluster_use = cluster_use,
    pw_param = pw_sample_stage,
    pathways = pathways,
    pairs_list = list_stage,
    save_dir = SAVE_DIR
  )
}

