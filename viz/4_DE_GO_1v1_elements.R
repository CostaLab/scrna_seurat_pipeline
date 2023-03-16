DEGO_1v1_elements <- function(scrna){

  cluster_use <- cluster
  dego_sample_name <- paste0("dego_name_", cluster_use)

  if(dego_sample_name %ni% names(scrna@tools)){
    stop(glue("ERROR: DE&GO samples comparing hasn't been calculated for cluster:{cluster_use}\n Please run [scrna_dego_name]!!!"))
  }



  dego_df <- seutools_partition(scrna,
                                partition=dego_sample_name,
                                save_dir=SAVE_DIR,
                                allinone=ALLINONE)
  all_de_list      <- dego_df$de
  all_goup_list    <- dego_df$goup
  all_godown_list  <- dego_df$godown

  sample_names <- scrna@tools$meta_order$name
  list_1v1 <- comb_list(sample_names)

  compare_vs_loop(
    cluster_use = cluster_use,
    all_de_list = all_de_list,
    all_goup_list = all_goup_list,
    all_godown_list = all_godown_list,
    pairs_list = list_1v1
  )
}

