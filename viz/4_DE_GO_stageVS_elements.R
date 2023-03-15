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

