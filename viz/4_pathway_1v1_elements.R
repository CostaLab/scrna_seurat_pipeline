pathway_1v1_elements <- function(scrna){
  cluster_use <- cluster

  dic <- c(
    "hallmark_1v1 " = "hallmark",
    "reactome_1v1" = "reactome",
    "kegg_1v1" = "kegg"
  )

  pathways <- dic[EXEC_PLAN]
  pathways <-pathways[!is.na(pathways)]

  pw_sample_name <- paste0("pathway_name_", cluster_use)

  if(pw_sample_name %ni% names(scrna@tools)){
    stop(glue("ERROR:{pw} samples comparing hasn't been calculated for cluster:{cluster_use}\n
              Please run [scrna_pathway_name]!!!"))
  }

  sample_names <- scrna@tools$meta_order$name
  list_1v1 <- comb_list(sample_names)

  pathway_vs_loop(
    scrna = scrna,
    cluster_use = cluster_use,
    pw_param = pw_sample_name,
    pathways = pathways,
    pairs_list = list_1v1
  )
}

