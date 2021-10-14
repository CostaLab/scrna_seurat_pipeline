DEGO_1v1_elements <- function(scrna){

  cluster_use <- cluster
  dego_sample_name <- paste0("dego_name_", cluster_use)

  if(dego_sample_name %ni% names(scrna@tools)){
    stop(glue("ERROR: DE&GO samples comparing hasn't been calculated for cluster:{cluster_use}\n Please run [scrna_dego_name]!!!"))
  }

  de.list <- scrna@tools[[dego_sample_name]]
  all_de_list <- scrna@tools[[dego_sample_name]]$de
  all_goup_list <- scrna@tools[[dego_sample_name]]$goup
  all_godown_list <- scrna@tools[[dego_sample_name]]$godown

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

