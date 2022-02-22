suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Libra))
suppressPackageStartupMessages(library(Seurat))

pseudobulk_de <- function(
  seurat_object,
  is_single_comparison = TRUE,
  replicate_col = "name",
  cell_type_col = "seurat_clusters",
  label_col = "stage",
  group_comparison_df = NULL
  ){

  n_groups <- seurat_object@meta.data[, label_col] %>% unique %>% length
  assertthat::assert_that(n_groups > 1)

  l_res <- NULL

  if(is_single_comparison){
    # 1 v 1 comparison
    # here we want to subset the obj to only compare group X against group Y
    groups_to_process <- group_comparison_df
    if(is.null(group_comparison_df)){
      groups_to_process <-
        seurat_object@meta.data %>%
        # make df with group combinations
        tidyr::expand(
          g1 = !!sym(label_col),
          g2 = !!sym(label_col)
          ) %>%
        # filtering wont work if g1 or g2 are factors
        dplyr::mutate(g1 = as.character(g1), g2 = as.character(g2)) %>%
        # filter out non-unique combinations
        dplyr::filter(g1 < g2)
    }

    # make sure groups defined are in fact in the data
    group1_in_data <- intersect(
      as.data.frame(groups_to_process)[,1] %>% unique,
      seurat_object@meta.data[, label_col] %>% unique
      ) %>% length
    group2_in_data <- intersect(
      as.data.frame(groups_to_process)[,2] %>% unique,
      seurat_object@meta.data[, label_col] %>% unique
      ) %>% length
    assertthat::assert_that(group1_in_data > 0)
    assertthat::assert_that(group2_in_data > 0)

    l_res <- purrr::map2(
      groups_to_process[,1],
      groups_to_process[,2],
      function(c1, c2){

        keep_samples <- 
          seurat_object@meta.data[, label_col] %in% c1 |
          seurat_object@meta.data[, label_col] %in% c2
        tmp_seurat <- seurat_object[, keep_samples]

        tmp_seurat@meta.data$tmp_label <-
          factor(tmp_seurat@meta.data[, label_col], levels = c(c1, c2))

        pseudo_de <- tmp_seurat %>%
          Libra::run_de(
            de_family = "pseudobulk",
            de_method = "limma",
            replicate_col = replicate_col,
            cell_type_col = cell_type_col,
            label_col = "tmp_label"
          )

        pseudo_de <- pseudo_de %>%
          dplyr::arrange(dplyr::desc(avg_logFC)) %>%
          dplyr::mutate(comparison = paste0(c1, "_", c2))
        return(pseudo_de)
    })

    names(l_res) <- paste0(groups_to_process[,1],"_v_",groups_to_process[,2])

  }else{
    # 1 v all comparison
    # here we want to effectively modify the label such that we can compare
    # group X against group A+B+C+Y+...
    groups_to_process <- 
      seurat_object@meta.data[, label_col] %>% unique %>% sort

    l_res <- purrr::map(groups_to_process,function(gg){

      tmp_label <- ifelse(seurat_object@meta.data[, label_col] == gg, gg, "others")

      seurat_object@meta.data$tmp_label <- factor(tmp_label, levels = c(gg, "others"))

      pseudo_de <- seurat_object %>%
        Libra::run_de(
          de_family = "pseudobulk",
          de_method = "limma",
          replicate_col = replicate_col,
          cell_type_col = cell_type_col,
          label_col = "tmp_label"
        )

      pseudo_de <- pseudo_de %>%
        dplyr::arrange(dplyr::desc(avg_logFC)) %>%
        dplyr::mutate(comparison = paste0(gg, "_", "others"))
      return(pseudo_de)
    })
    names(l_res) <- paste0(groups_to_process, "_v_", "others")
  }

  return(l_res)
}


# Examples ####
if(FALSE){

  seurat_obj <- CimpleG::load_object("results/exampleproj/save/scrna_phase_comparing.Rds")

  seurat_obj$label <-
    sample(
      c("AA", "BB", "CC", "DD"),
      seurat_obj@meta.data%>%nrow(),
      replace=TRUE,
      prob = c(.5,.3,.15,.15)
    )

  # relevant for scDE methods
  calc_dv <- Libra::calculate_delta_variance(
    seurat_obj,
    replicate_col = "name",
    cell_type_col = "seurat_clusters",
    label_col="label"
  )

  ### name -------------- replicate
  ### seurat_clusters --- cell_type
  ### stage ------------- label

  # we want to make a function to compare whichever groups of 2 in a vector
  # Cond1 vs Cond2 / Cond1 vs Cond3 / etc
  # we want to make a function to compare 1 group vs all others
  # Cluster1 vs ALL / Cluster2 vs All / etc

  example_1v1 <- pseudobulk_de(
    seurat_obj,
    replicate_col = "name",
    cell_type_col = "seurat_clusters",
    label_col="stage",
    is_single_comparison=TRUE
  )

  example_1v1_explicitgroup.1 <- pseudobulk_de(
    seurat_obj,
    replicate_col = "name",
    cell_type_col = "seurat_clusters",
    label_col="stage",
    is_single_comparison=TRUE,
    group_comparison_df = data.frame(g1="MxCre",g2="Csnk")
  )

  example_1v1_explicitgroup.2 <- pseudobulk_de(
    seurat_obj,
    replicate_col = "name",
    cell_type_col = "seurat_clusters",
    label_col="label",
    is_single_comparison=TRUE,
    group_comparison_df = data.frame(g1=c("CC","AA"),g2=c("DD","CC"))
  )


  seurat_obj@meta.data$cell_type <- "bulk"

  example_1vAll <- pseudobulk_de(
    seurat_obj,
    replicate_col = "name",
    cell_type_col = "cell_type",
    label_col="seurat_clusters",
    is_single_comparison=FALSE
  )


  get_top_n <- 10
  example_1v1[[1]] %>% dplyr::pull(cell_type) %>% table

  example_1v1[[1]] %>%
    dplyr::filter(cell_type == "0") %>%
    #     dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(dplyr::desc(avg_logFC)) %>%
    dplyr::select(cell_type,gene,avg_logFC,p_val_adj,comparison) %>%
    #     dplyr::filter(gene=="mt-Atp6")
    dplyr::slice(1:get_top_n,(dplyr::n()-(get_top_n-1)):dplyr::n())

  example_1vAll[[4]] %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(dplyr::desc(avg_logFC)) %>%
    dplyr::slice(1:get_top_n,(dplyr::n()-(get_top_n-1)):dplyr::n())

  str(asd)
}


