#!/usr/bin/Rscript

t_start <- Sys.time()

source("R/helper_functions.R")
source("R/DE_GO_VS_helper.R")
source("R/pathway_vs_helper.R")
source("R/save_load_helper.R")
source("R/scProportion.R")


suppressPackageStartupMessages(library(optparse))      ## Options

###======================PARAMETERS BEGIN====================================
AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(
    parser, c("-a", "--author"), type = "character", default = "Costalab",
    help = "Author display in report [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-m", "--make_element"), type = "character", default = "FALSE",
    help = "make_element [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-s", "--savedir"), type = "character", default = "save",
    help = "savedir [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-c", "--configfile"), type = "character", default = "conf/config.R",
    help = "configfile [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-o", "--report_dir"), type = "character", default = "report",
    help = "report_dir [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-r", "--charts_dir"), type = "character", default = "charts",
    help = "charts_dir [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-p", "--project"), type = "character", default = "",
    help = "charts_dir [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-f", "--featureplotstyle"), type = "character", default = "seurat",
    help = "feature style: seurat, schex, nebulosa, current  [default %default]", metavar = "character")

  parser <- add_option(
    parser, c("-e", "--externalfile"), type = "character",
    default = "./external/Human_and_mouse_cell_markers-Markers.tsv",
    help = "report_dir [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-d", "--defaultclsuters"), type = "character",
    default = "seurat_clusters",
    help = paste0(
      "cluster slots: \n",
      "removed_clusters\n",
      "remove_reclusters\n",
      "merged_clusters\n",
      "annotation\n",
      "singleton [default %default]"),
    metavar = "character")
  parser <- add_option(
    parser, c("-j", "--planOfreport"), type = "character",
    default = "[\"DEGO_1v1\"]",
    help = "plan of report [default %default]", metavar = "character")
  parser <- add_option(
    parser, c("-l", "--singlefile"), type = "character", default = "FALSE",
    help = "generate single html file [default %default]", metavar = "character")

  parser <- add_option(
    parser, c("-i", "--indexonly"), type = "character", default = "FALSE",
    help = "only generate index.html [default %default]", metavar = "character")

  parser <- add_option(
    parser, c("-z", "--compression"),
    type = "character", default = "gzip", metavar = "character",
    help = paste0(
      "Compression algorithm to use when saving R objects [default %default]. ",
      "One of 'zstd', 'lz4', 'gzip', 'bzip2', ",
      "'xz' or 'nocomp' (no compression)"
    )
  )
  return(parser)
}

print_nelement_msg <- function(element){
  cat(
    paste0(
      date(),
      green(" Making element: "),
      red(as.character(element)), "\n"
    )
  )
}
print_ngenelement_msg <- function(element){
  cat(
    paste0(
      date(),
      blue(" Generating: "),
      red(as.character(element)), "\n"
    )
  )
}

parser <- AllOptions()
pa <- parse_args(parser)

PROJECT           = pa$project
AUTHOR            = pa$author
MAKE_ELEMENT      = pa$make_element
SAVE_DIR          = pa$savedir
CONFIGFILE        = pa$configfile
REPORTDIR         = pa$report_dir
CHARTSDIR         = pa$charts_dir
EXTERNALFILE      = pa$externalfile
DEFAULTCLUSTERS   = pa$defaultclsuters
GEN_SINGLE_FILE   = pa$singlefile
EXEC_PLAN         = jsonlite::fromJSON(pa$planOfreport)
INDEX_ONLY        = pa$indexonly
COMPRESSION_FORMAT= pa$compression
FEATUREPLOT_STYLE = pa$featureplotstyle
dir.create(file.path(REPORTDIR, "data"), recursive = TRUE)


py_exe_list = paste0(EXEC_PLAN, collapse = "','")
py_exe_list = paste0("['", py_exe_list, "']")
###=======================PARAMETERS END=============================

message("Project:", PROJECT)
cat(paste("Use cluster slot", red(DEFAULTCLUSTERS), "\n"))
##1. copy excels to report/data
message("copying excel files from ", CHARTSDIR, "\n")
file.copy(
  file.path(CHARTSDIR, list.files(CHARTSDIR)),
  glue("{REPORTDIR}/data"),
  overwrite = TRUE, copy.date = TRUE
)
dir.create(glue("{REPORTDIR}/../viz"))
file.copy(
  file.path(getwd(), "viz", list.files("viz")),
  glue("{REPORTDIR}/../viz"),
  overwrite = TRUE, recursive = TRUE, copy.date = TRUE
)

code_gene_file <- file.path(REPORTDIR, "..", "viz/code_generator.py")

viz_path <- file.path(REPORTDIR, "..", "viz")

code_generate_cmd <- glue(
  "python {code_gene_file} ",
  "-c {DEFAULTCLUSTERS} ",
  "-cf '{CONFIGFILE}' ",
  "--save_dir '{SAVE_DIR}' ",
  "--output_dir '{REPORTDIR}' ",
  "--proj_tag {PROJECT} ",
  "--executing_list \"{py_exe_list}\""
)


##2. code generate
run_shell(code_generate_cmd)
##3. generate index.html from template
run_shell(glue("grip --export {REPORTDIR}/index.md"))

if(INDEX_ONLY){
  cat(red("======Only generate index.html=====\n"))
  quit(save = "no")
}

report_data_folder  = file.path(REPORTDIR, "data")
report_tables_folder = file.path(REPORTDIR, "tables")
report_plots_folder = file.path(REPORTDIR, "plots")
report_plots_folder_png = file.path(report_plots_folder, "png")
report_plots_folder_pdf = file.path(report_plots_folder, "pdf")
savedir = SAVE_DIR
cluster = DEFAULTCLUSTERS
project = PROJECT
funcs = EXEC_PLAN
ext_annot_fp = EXTERNALFILE


##4. Make_report element
source(CONFIGFILE)
scrna <- NULL


if(MAKE_ELEMENT){

  if(identical(cluster,"singleton")){
    cat(
      paste0(date(), blue(" Loading: "), red("scrna_phase_singleton.Rds"), "\n")
    )
    scrna <- load_object(
      file_name = file.path(savedir, "scrna_phase_singleton.Rds")
    )
    cat(
      paste0(date(), blue(" Loaded: "), red("scrna_phase_singleton.Rds"), "\n")
    )
  }else{
    cat(
      paste0(date(), blue(" Loading: "), red("scrna_phase_comparing.Rds"), "\n")
    )
    scrna <- load_object(
      file_name = file.path(savedir, "scrna_phase_comparing.Rds")
    )
    cat(
      paste0(date(), blue(" Loaded: "), red("scrna_phase_comparing.Rds"), "\n")
    )
  }

  scrna$name <- factor(scrna$name, levels = names(data_src))
  scrna$stage <- factor(scrna$stage, levels = unique(stage_lst))
  cluster_viridis_opt = ifelse(
    any(grepl("cluster_color_option", names(viz_conf), fixed = TRUE)),
    viz_conf[["cluster_color_option"]], # Config option
    "d3" # Default
  )

  # replicate_colors
  replicates_viridis_opt = ifelse(
    any(grepl("replicate_color_option", names(viz_conf), fixed = TRUE)),
    viz_conf[["replicate_color_option"]],
    "futurama"
  )
  neg_color = ifelse(
    any(grepl("neg_color", names(viz_conf), fixed = TRUE)),
    viz_conf[["neg_color"]],
    "#51C3CC" # colorBlindness::Blue2DarkOrange12Steps[2]
  )

  # divergent color for positive values
  pos_color = ifelse(
    any(grepl("pos_color", names(viz_conf), fixed = TRUE)),
    viz_conf[["pos_color"]],
    "#CC5800" #rev(colorBlindness::Blue2DarkOrange12Steps)[2]
  )

  # divergent color for base values
  base_color = ifelse(
    any(grepl("base_color", names(viz_conf), fixed = TRUE)),
    viz_conf[["base_color"]],
    "lightgrey"
  )

  # divergent palette for neg to pos values
  neg_pos_divergent_palette =
    if(any(grepl("neg_pos_divergent_palette", names(viz_conf), fixed = TRUE))){
      viz_conf[["neg_pos_divergent_palette"]]
    }else{
      # colorBlindness::Blue2DarkOrange12Steps
      c("#1E8E99", "#51C3CC", "#99F9FF", "#B2FCFF",
        "#CCFEFF", "#E5FFFF", "#FFE5CC", "#FFCA99",
        "#FFAD65", "#FF8E32", "#CC5800", "#993F00")
    }

  zero_pos_divergent_colors = c(base_color, pos_color)

  override_colours = c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
    "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99"
  )



  available_clusters = DEFAULTCLUSTERS
  if("Clusters_harmony" %in% EXEC_PLAN){
    available_clusters <- c(available_clusters, "harmony_inte_clusters")
  }
  if("Clusters_seurat" %in% EXEC_PLAN){
    available_clusters <- c(available_clusters, "seurat_inte_clusters")
  }


  dir.create(report_plots_folder_png, recursive = TRUE)
  dir.create(report_plots_folder_pdf, recursive = TRUE)
  dir.create(report_tables_folder, recursive = TRUE)


  source(glue("{viz_path}/1_quality_report_elements.R"))
  source(glue("{viz_path}/ambientRNA_viz_elements.R"))
  source(glue("{viz_path}/doubletdetection_viz_elements.R"))
  source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  source(glue("{viz_path}/2_batch_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clustering_elements.R"))
  source(glue("{viz_path}/2_clusters_DEs_elements.R"))
  source(glue("{viz_path}/3_external_markers_elements.R"))
  source(glue("{viz_path}/3_DE_GO-analysis_elements.R"))
  source(glue("{viz_path}/4_DE_GO_1v1_elements.R"))
  source(glue("{viz_path}/4_DE_GO_stageVS_elements.R"))
  source(glue("{viz_path}/4_pathway_1v1_elements.R"))
  source(glue("{viz_path}/4_pathway_stageVS_elements.R"))
  source(glue("{viz_path}/4_Genesets_1v1_elements.R"))
  source(glue("{viz_path}/4_Genesets_stageVS_elements.R"))
  source(glue("{viz_path}/4_progeny_stageVS_elements.R"))




# run necessary generators
  if("QC" %in% EXEC_PLAN) {
    print_nelement_msg("QC")
    quality_report_elements() ## load scrna innner the function
  }


  if(any(grepl("AmbientRNA", funcs, fixed = TRUE))){
    print_nelement_msg("AmbientRNA")
    ambientRNA_elements(scrna)
  }
  if(any(grepl("DoubletDetection",funcs,fixed=TRUE))){
    print_nelement_msg("DoubletDetection")
    doubletdetection_viz_elements(scrna)
  }
  if("DEs"%in% EXEC_PLAN) {
    print_nelement_msg("DEs")
    clusters_DEs_elements(scrna)
  }
  if(any(grepl("Clusters_", EXEC_PLAN, fixed = TRUE))){
    print_nelement_msg("Clusters_integration")
    batch_clustering_elements(scrna)
  }
  if("Clusters" %in% EXEC_PLAN){
    print_nelement_msg("Clusters")
    clustering_elements(scrna)
  }
  if("Singleton" %in% EXEC_PLAN){
    print_nelement_msg("Singleton")
    clustering_elements(scrna)
    clusters_DEs_elements(scrna)
  }
  if("EXT_MARKERS" %in% EXEC_PLAN){
    print_nelement_msg("EXT_MARKERS")
    external_markers_elements(scrna)
  }
  # FIXME possible problem where all term enrichment analysis is on the same place
  if(length(intersect(c("DEGO","Genesets","progeny","hallmark","KEGG","Reactome"), EXEC_PLAN) > 0)){
    print_nelement_msg("DEGO")
    DE_GO_analysis_elements(scrna)
  }
  if("DEGO_1v1" %in% EXEC_PLAN){
    print_nelement_msg("DEGO_1v1")
    DEGO_1v1_elements(scrna)
  }
  if("DEGO_stage" %in% EXEC_PLAN){
    print_nelement_msg("DEGO_stage")
    DEGO_stageVS_elements(scrna)
  }
  if(length(intersect(c("hallmark_1v1","reactome_1v1","kegg_1v1"), EXEC_PLAN) > 0)){
    print_nelement_msg("pathway_1v1")
    pathway_1v1_elements(scrna)
  }
  if(length(intersect(c("hallmark_stage","reactome_stage","kegg_stage"), EXEC_PLAN) > 0)){
    print_nelement_msg("pathway_stage")
    pathway_stage_elements(scrna)
  }
  if(length(intersect(c("Genesets_stage"), EXEC_PLAN) > 0)){
    print_nelement_msg("Genesets_stage")
    Genesets_stageVS_elements(scrna)
  }
  if(length(intersect(c("Genesets_1v1"), EXEC_PLAN) > 0)){
    print_nelement_msg("Genesets_1v1")
    Genesets_1v1_elements(scrna)
  }
  if(length(intersect(c("progeny_stage"), EXEC_PLAN) > 0)){
    print_nelement_msg("progeny_stage")
    progeny_stageVS_elements(scrna)
  }
}

cluster_info <- build_cluster_info(scrna)

##5. Produce Report
render_func = function(rmd_input_filename, output_filename){
  rmarkdown::render(
    rmd_input_filename,
    output_file = output_filename,
    output_dir = file.path(REPORTDIR,"data"),
    output_format = c("html_document"),
    quiet = TRUE,
    clean = TRUE,
    params = list(
      scrna = scrna,
      cluster = DEFAULTCLUSTERS,
      cluster_info = cluster_info,
      project = PROJECT,
      savedir = SAVE_DIR,
      funcs = EXEC_PLAN,
      report_data_folder = report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author = AUTHOR
    )
  )
}

dic_Rmd_n_Output <- list(
  "QC"               = c(glue("{viz_path}/1_quality_report.Rmd"),        "data_quality"),
  "AmbientRNA"       = c(glue("{viz_path}/ambientRNA_viz.Rmd"),          "ambient_rna"),
  "DoubletDetection" = c(glue("{viz_path}/doubletdetection_viz.Rmd"),    "doublet_detection"),
  "DEs"              = c(glue("{viz_path}/2_clusters_DEs.Rmd"),          "clusters_DEs"),
  "Clusters"         = c(glue("{viz_path}/2_clustering.Rmd"),            "clusters"),
  "Singleton"        = c(glue("{viz_path}/2_clustering.Rmd"),            "clusters"),
  "Clusters_harmony" = c(glue("{viz_path}/2_clustering_harmony.Rmd"),    "clusters_harmony"),
  "Clusters_seurat"  = c(glue("{viz_path}/2_clustering_seurat.Rmd"),     "clusters_seurat"),
  "EXT_MARKERS"      = c(glue("{viz_path}/3_external_markers.Rmd"),      "external_markers"),
  "DEGO"             = c(glue("{viz_path}/3_DE_GO-analysis.Rmd"),        "dego"),
  "KEGG"             = c(glue("{viz_path}/3_KEGG.Rmd"),                  "KEGG"),
  "progeny"          = c(glue("{viz_path}/3_progeny.Rmd"),               "progeny"),
  "Genesets"         = c(glue("{viz_path}/3_Genesets.Rmd"),              "Genesets"),
  "hallmark"         = c(glue("{viz_path}/3_hallmark.Rmd"),              "hallmark"),
  "Reactome"         = c(glue("{viz_path}/3_Reactome.Rmd"),              "Reactome"),
  "DEGO_stage"       = c(glue("{viz_path}/4_DE_GO_%s.vs.%s_stageVS.Rmd"),"gv"),
  "DEGO_1v1"         = c(glue("{viz_path}/4_DE_GO_%s.vs.%s_1v1.Rmd"),    "1vs1"),
  "hallmark_1v1"     = c(glue("{viz_path}/4_hallmark_1v1.Rmd"),          "hallmark_1vs1"),
  "reactome_1v1"     = c(glue("{viz_path}/4_reactome_1v1.Rmd"),          "reactome_1vs1"),
  "kegg_1v1"         = c(glue("{viz_path}/4_kegg_1v1.Rmd"),              "kegg_1vs1"),
  "hallmark_stage"   = c(glue("{viz_path}/4_hallmark_stageVS.Rmd"),      "hallmark_stageVS"),
  "reactome_stage"   = c(glue("{viz_path}/4_reactome_stageVS.Rmd"),      "reactome_stageVS"),
  "kegg_stage"       = c(glue("{viz_path}/4_kegg_stageVS.Rmd"),          "kegg_stageVS"),
  "Genesets_1v1"     = c(glue("{viz_path}/4_Genesets_1v1.Rmd"),          "Genesets_1vs1"),
  "Genesets_stage"   = c(glue("{viz_path}/4_Genesets_stageVS.Rmd"),      "Genesets_stageVS"),
  "progeny_stage"    = c(glue("{viz_path}/4_progeny_stageVS.Rmd"),       "progeny_stageVS"),
  "intUMAPs"         = c(glue("{viz_path}/interactive_UMAPs.Rmd"),       "interactive_UMAPs")
)




for(exec_elem in EXEC_PLAN){
  print_ngenelement_msg(exec_elem)
  if(exec_elem %ni% names(dic_Rmd_n_Output)){
    message("viz ", exec_elem, " is not implemented!!")
    next
  }
  rmd_n_output <- dic_Rmd_n_Output[[exec_elem]]
  rmd <- rmd_n_output[1]
  output <- rmd_n_output[2]
  if(output == "1vs1"){
    for(apair in comb_list(names(data_src))){
      render_func(
        sprintf(rmd, apair[1], apair[2]),
        glue("{output}_{apair[1]}.vs.{apair[2]}.html")
      )
    }
  }else if(output == "gv"){
    for(apair in comb_list(unique(stage_lst))){
      render_func(
        sprintf(rmd, apair[1], apair[2]),
        glue("{output}_{apair[1]}.vs.{apair[2]}.html")
      )
    }
  }else{
    render_func(rmd, output)
  }

}

if(GEN_SINGLE_FILE){
  # generate report including everything in a single file
  rmarkdown::render(
    glue("{viz_path}/scrna_pipeline_report.Rmd"),
    output_file = file.path(REPORTDIR, "scrna_report"),
    output_format = c("html_document"),
    clean = TRUE,
    params = list(
      cluster = DEFAULTCLUSTERS,
      project = PROJECT,
      savedir = SAVE_DIR,
      funcs = EXEC_PLAN,
      report_data_folder = report_data_folder,
      report_tables_folder = report_tables_folder,
      report_plots_folder = report_plots_folder,
      report_plots_folder_png = report_plots_folder_png,
      report_plots_folder_pdf = report_plots_folder_pdf,
      author = AUTHOR
    )
   )
}

t_diff <- Sys.time() - t_start
cat(
  black$bgGreen$bold(paste0(
      "Time difference of ", as.character(round(t_diff, 2)), " ", units(t_diff)
  )), "\n"
)

