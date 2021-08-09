DEGO_stageVS_elements <- function(scrna){

  cluster_use <- cluster
  dego_stage_name <- paste0("dego_stage_", cluster_use)

  if(dego_stage_name %ni% names(scrna@tools)){
    stop(glue("ERROR: DE&GO stages comparing hasn't been calculated for cluster:{cluster_use}\n Please run [scrna_dego_name]!!!"))
  }

  de.list <- scrna@tools[[dego_stage_name]]
  all_de_list <- scrna@tools[[dego_stage_name]]$de
  all_goup_list <- scrna@tools[[dego_stage_name]]$goup
  all_godown_list <- scrna@tools[[dego_stage_name]]$godown

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

