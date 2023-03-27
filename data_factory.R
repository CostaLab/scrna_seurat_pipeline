#!/usr/bin/Rscript

###Set VERSION
VERSION = "1.0.5"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
  message("scRNA seurat pipeline\nVersion: \n\t", VERSION, "\n")
  quit(status = 0)
}


suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger
suppressPackageStartupMessages(library(dplyr))
`%ni%` <- Negate(`%in%`)




### TODO LIST
#--- save the parameter into R object
#--- 1v1 cluster name set as a parameter should be good

AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-n", "--numberofcores"), type="integer", default=1,
                       help="Parallel Run Number of cores [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--MaxMemMega"), type="integer", default=20000,
                       help="Parallel Max Memory size megabytes [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("-c", "--configfile"), type="character", default="conf/config.R",
                       help="Config file for the run, input files and executing plan settings  [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("-l", "--logname"), type="character", default="",
                       help="add a unique name for log files [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("--nFeatureRNAfloor"), type="character", default=400,
                       help="min nFeature_RNA [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--nFeatureRNAceiling"), type="character", default=Inf,
                       help="max nFeature_RNA [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--nCountRNAfloor"), type="character", default=0,
                       help="min nCount_RNA [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--nCountRNAceiling"), type="character", default=40000,
                       help="max nCount_RNA [default %default]",
                       metavar="number")
  parser <- add_option(parser, c("--pct_mitofloor"), type="character", default=-Inf,
                       help="min percentage mito [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--pct_mitoceiling"), type="character", default=100,
                       help="max percentage mito [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--pct_ribofloor"), type="character", default=-Inf,
                       help="min percentage ribo [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("--pct_riboceiling"), type="character", default=100,
                       help="max percentage ribo [default %default]",
                       metavar="number")

  parser <- add_option(parser, c("-f", "--countmatrixformat"), type="character", default="10X",
                       help="count matrix format, current support 10X and 10X_h5  [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("-s", "--savedir"), type="character", default="save",
                       help="RObject save directory [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("-e", "--chartsdir"), type="character", default="charts",
                       help="xlsx directory [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("-d", "--dims4Integrate"), type="character", default="1:30",
                       help="Dims to keep for integrating [default %default]",
                       metavar="VECTOR")

  parser <- add_option(parser, c("-x", "--Dims4FindNeighbors"), type="character", default="1:50",
                       help="Dims kept for findNeighbors [default %default]",
                       metavar="VECTOR")

  parser <- add_option(parser, c("-r", "--clusterresolution"), type="numeric", default=0.5,
                       help="Resolution for clustering [default %default]",
                       metavar="numeric")

  parser <- add_option(parser, c("-a", "--defaultclustername"), type="character", default="seurat_clusters",
                       help="downstream analysis default cluster name  [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("--harmony_dim"), type="character", default="1:50",
                       help="dims keep to do clustering using harmony [default %default]",
                       metavar="character")

  parser <- add_option(parser, c("--allinone"), type="logical", default=TRUE,
                       help="store calculated result in scrna@tools & assays in scrna@assays if TRUE, else store in save/partition & save/assays respectively [default %default]",
                       metavar="character")

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
parser <- AllOptions()
#debug(parse_args)
pa <- parse_args(parser)



WORKER_NUM            = pa$numberofcores
MAXMEMMEGA            = pa$MaxMemMega
SAVE_DIR              = pa$savedir
CHARTS_DIR            = pa$chartsdir
INTEGRATED_DIM        = eval(parse(text=pa$dims4Integrate))
HARMONY_DIM           = eval(parse(text=pa$harmony_dim))
FINDNEIGHBORS_DIM     = eval(parse(text=pa$Dims4FindNeighbors))
CLUSTER_RESOLUTION    = pa$clusterresolution
DEFUALT_CLUSTER_NAME  = pa$defaultclustername
CM_FORMAT             = pa$countmatrixformat
COMPRESSION_FORMAT    = pa$compression
ALLINONE              = pa$allinone


f = function(x){if(x == ""){ return("")}else{ return("-")}}
LOGNAME               = paste0(pa$logname, f(pa$logname))


toInt <- function(x) switch(x, "Inf"=Inf, "-Inf" = -Inf, as.integer(x))
nFeatureRNAfloor      = toInt(pa$nFeatureRNAfloor)
nFeatureRNAceiling      = toInt(pa$nFeatureRNAceiling)
nCountRNAceiling        = toInt(pa$nCountRNAceiling)
nCountRNAfloor       = toInt(pa$nCountRNAfloor)
pct_mitoceiling       = toInt(pa$pct_mitoceiling)
pct_mitofloor         = toInt(pa$pct_mitofloor)
pct_riboceiling       = toInt(pa$pct_riboceiling)
pct_ribofloor         = toInt(pa$pct_ribofloor)


CLUSTER_RESOLUTION_RANGE = seq(0.1, 0.8, 0.1)


if(!file.exists("logs")){
  dir.create("logs")
}

cur_date <- as.character(Sys.Date())
errorLog  =  file.path("logs", sprintf("%sERROR-%s.log", LOGNAME, cur_date))
warnLog   =  file.path("logs", sprintf("%sWARN-%s.log", LOGNAME, cur_date))
infoLog   =  file.path("logs",  sprintf("%sINFO-%s.log", LOGNAME, cur_date))
traceLog   =  file.path("logs",  sprintf("%sTRAC-%s.log", LOGNAME, cur_date))

invisible(flog.logger("error", ERROR, appender.file(errorLog)))
invisible(flog.logger("warn", WARN, appender.file(warnLog)))
invisible(flog.logger("info", INFO, appender.file(infoLog)))
invisible(flog.logger("trace", TRACE, appender.file(traceLog)))
invisible(flog.appender(appender.console(), name = "ROOT"))

logger.info <- function(msg, ...) {
  flog.info(msg, ..., name = "ROOT")
  flog.info(msg, ..., name = "info")
}

logger.warn <- function(msg, ...) {
  flog.warn(msg, ..., name = "ROOT")
  flog.warn(msg, ..., name = "info")
  flog.warn(msg, ..., name = "warn")
}

logger.error <- function(msg, ...) {
  flog.error(msg, ..., name = "ROOT")
  flog.error(msg, ..., name = "info")
  flog.error(msg, ..., name = "warn")
  flog.error(msg, ..., name = "error")
  flog.error(msg, ..., name = "trace")
}


## Robject directory
if(!file.exists(SAVE_DIR)){
  dir.create(SAVE_DIR)
  logger.info("RObject save in directory: %s", SAVE_DIR)
}else{
  logger.info("RObject save in existed directory: %s", SAVE_DIR)
  logger.info("meta.data save in directory: %s", file.path(SAVE_DIR))
}

if(!file.exists(file.path(SAVE_DIR, "meta"))){
  dir.create(file.path(SAVE_DIR, "meta"), showWarnings=F)
  logger.info("meta.data save in directory: %s", file.path(SAVE_DIR, "meta"))
}

if(!ALLINONE){
    dir.create(file.path(SAVE_DIR, "partition"), showWarnings=F)
    logger.info("data in @tools save in directory: %s/partition", SAVE_DIR)
    dir.create(file.path(SAVE_DIR, "assays"), showWarnings=F)
    logger.info("data in assays save in directory: %s/assays", SAVE_DIR)
}


if(!file.exists(CHARTS_DIR)){
  dir.create(CHARTS_DIR)
  logger.info("xlsx save in directory: %s", CHARTS_DIR)
}else{
  logger.info("xlsx save in existed directory: %s", CHARTS_DIR)
}


## configuration file
if(!file.exists(pa$configfile)){
  logger.error("Please check the config file!")
  print("Please check the config file!")
  parse_args(parser, args = c("--help"))
}

source(pa$configfile)
source("R/save_load_helper.R")

##--------------Load the function scripts----------------------
R_scripts <- list.files("R_scripts/", full.name = TRUE)
for(i in 1:length(R_scripts)) source(R_scripts[i])

##--------------keep variables----------------------
pa$project_name       = PROJECT
pa$organ              = ORGAN
pa$species            = SPECIES
pa$mca_name           = MCA_NAME
pa$hcl_name           = HCL_NAME
pa$external_file      = ANNOTATION_EXTERNAL_FILE
pa$data_src           = data_src
pa$stage_lst          = stage_lst
pa$integration_option = INTEGRATION_OPTION
pa$doublet_switch     = doublet_switch
pa$doublet_lst        = doublet_lst

##--------------Loading the relevant slots and attributes for the sanity function----------------------
sanity_attributes <- read.table("static/SlotsAttributes.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)


##--------------Doublet estimation----------------------
estimation_10X <- read.csv("static/DoubletEstimation10X.csv")
model_recovered <- approxfun(x = estimation_10X$CellsRecovered, y = estimation_10X$MultipletRate, rule = 2)


logger.info("----------------------Project: %s --------------------------", PROJECT)


## src data summary
logger.info("samples summary: ")
for(i in seq_along(data_src)){
  key = names(data_src)[i]
  val = data_src[[key]]
  logger.info(paste(i, key, val, sep= "   ------->  "))
}



plan_record <- function(plan){
  start_line <- "\n\n======================================================================="
  rfile <- "logs/plans.log"
  write(start_line, rfile, append=TRUE)
  write(as.character(Sys.time()), rfile, append=TRUE)
  write("Parameters:", rfile, append=TRUE)
  for (nm in names(pa)) {
    write(sprintf("\t\t\t%s:  %s", nm, pa[nm]), rfile, append=TRUE)
  }
  write("\nplan:", rfile, append=TRUE)

  for (x in plan) {
    write(x, rfile, append=TRUE)
  }
}

# setup function to save object to debug for ease of use
debug_save_stop <- function(scrna,key){
  save_object(
    object = scrna,
    file_name = file.path(SAVE_DIR, "scrna_for_debug.Rds"),
    file_format = COMPRESSION_FORMAT
  )
  stop(glue("ERROR when run {key}"))
}

## executing plan
conf = conf[conf > 0]

if(conf[1] != 2 & (names(conf)[1] %ni% c("scrna_phase_preprocess", "scrna_phase_singleton", "scrna_rawdata"))){
  logger.error("!!Run without loading data, please check execution plan")
  stop("exit 1")
}

orig_len <- length(conf)
remove_useless_loadings <- function(conf){
  leng <- length(conf)
  if(conf[1] != 2){ ## no loading
    return(conf[1:leng])
  }

  for(i in 1:length(conf)){
    if(conf[i] == 2 ){  next  }
    else{return(conf[(i-1):leng]) }
  }
  return(conf[1:leng])
}

conf <- remove_useless_loadings(conf)
if(orig_len > length(conf)){
  logger.warn("Conf set too many loading Rds, will not load useless ones!")
}


dic = c("1"="calculating", "2"="loading")
plan = sprintf("\t\t\t%s: %s", dic[as.character(conf)], names(conf))
logger.info("The executing plan is: ")
logger.info(plan)
plan_record(plan)


## TODO program...
## Just One Sample should be taken into account
# library roxygen2 -> to maintain documents

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(progeny))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(future.apply))    ## parallel
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(celda))
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Rmagic))
suppressPackageStartupMessages(library(DoubletFinder))

registerDoParallel(cores=WORKER_NUM)


# Set up future for parallelization
plan("multicore", workers = WORKER_NUM)
options(future.globals.maxSize = MAXMEMMEGA * 1024^2)



cur_time <- as.character(Sys.time())
conf_main <- function(){
  scrna = NULL
  cluster.de = NULL
  for(i in seq_along(conf)){
    key = names(conf)[i]
    val = conf[i]
    NRds = paste(key, ".Rds", sep="")

    if(val==0){
      # do nothing
    }else if(val==2){
      logger.info(paste("loading", NRds))
      scrna <- load_object(file_name = file.path(SAVE_DIR, NRds))
      scrna@tools$parameter[[cur_time]] <- unlist(pa)
      scrna@tools$execution[[cur_time]] <- conf
      scrna <- sanity_function(scrna, key)
      scrna@tools$allinone = ALLINONE
      logger.info(paste("finished loading", NRds))
    }else if(val==1){
      func_name = paste("generate_", key, sep="")
      func_call = paste(func_name, "(scrna)", sep="")
      logger.info(paste("executing", func_name))
      if(startsWith(key, "scrna")){
        ret_list <-  eval(parse(text=func_call))
        scrna  <- ret_list[[1]]
        ret_code <- ret_list[[2]]
        scrna <- sanity_function(scrna, key)

        if(ret_code != 0){ debug_save_stop(scrna, key) }

        #if("meta_order" %ni% names(scrna@tools)){
        scrna@tools[["meta_order"]] <- list(
          name = names(data_src),
          stage = unique(stage_lst)
        )
        scrna@tools$allinone = ALLINONE
        #}
        phase_name <- get_output_name(scrna, key) ### only first store, or second time
        save_object(
          object = scrna,
          file_name = file.path(SAVE_DIR, glue("{phase_name}.Rds")),
          file_format = COMPRESSION_FORMAT
        )
        save_object(
          object = scrna@meta.data,
          file_name = file.path(SAVE_DIR, 'meta', glue("{phase_name}_meta.Rds")),
          file_format = COMPRESSION_FORMAT
        )


        logger.info(paste("finished", phase_name))
      }
    }else{
      logger.error("Wrong settings in config file in section RUN PARAMETERS\n Only 0, 1 or 2 permitted")
      logger.error(traceback())
      stop("Exit 1")
    }
  }
  logger.info("===============Finished===============")
}


get_regressout_vector <- function(){
    dic <- list(
            "mito"       = "percent.mt",
            "ribo"       = "percent.ribo",
            "cellcycle"  = c("G2M.Score", "S.Score"))

    keep <- names(preprocess_regressout[preprocess_regressout==1])
    ro <- unlist(dic[keep])
    names(ro) <- NULL
    return(ro)
}

get_output_name<- function(scrna, scrna_run_name){
    ret_code <- 0
    pconf <- configr::read.config("static/phase.ini")
    phases <- names(pconf)

    if(scrna_run_name %in% paste0("scrna_", phases)){
      return(scrna_run_name)
    }else if(length(unique(scrna$name)) == 1){
      return("scrna_phase_singleton")
    }else{
      return("scrna_phase_comparing")
    }

}


dego_dump <- function(file_predix, de.list, go_ups, go_downs) {
  #file_predix <- paste0(a_pair[1], ".vs.", a_pair[2])

  flist <- lapply(de.list, subset, subset = p_val_adj < 0.05)
  if(length(flist) == 0){
    return()
  }
  flist <-  flist[sapply(flist, function(m) nrow(m) >0)]
  WriteXLS(
           flist,
           file.path(CHARTS_DIR,  sprintf("%s.de_%s.xlsx", file_predix, DEFUALT_CLUSTER_NAME)),
           SheetNames = reformat_sheetnames(names(flist)))

  golist_xls(go_ups, sprintf("%s.goup_%s.xlsx",file_predix, DEFUALT_CLUSTER_NAME))
  golist_xls(go_downs,sprintf("%s.godown_%s.xlsx",file_predix, DEFUALT_CLUSTER_NAME))
}

pathway_dump <- function(file_predix, ptype, ups, downs) {
  golist_xls(ups, glue("{file_predix}.{ptype}up_{DEFUALT_CLUSTER_NAME}.xlsx"))
  golist_xls(downs,glue("{file_predix}.{ptype}down_{DEFUALT_CLUSTER_NAME}.xlsx"))
}

list_1vs1 <- function(){
  m <- combn(names(data_src), 2)
  n <- length(m)/2
  lst <- vector("list", n)

  for (i in 1:n){
    lst[[i]] <- m[1:2, i]
  }
  return(lst)
}


list_stages <- function(){
  conds <- unique(stage_lst)
  m <- combn(conds, 2)
  n = length(m)/2
  lst <- vector("list", n)
  for (i in 1:n){
    lst[[i]] <- m[1:2, i]
  }
  return(lst)
}


golist_xls <- function(go.list, fxls){
  suppressPackageStartupMessages(require(ggplot2))
  lst = lapply(go.list, function(x){
                 if(length(x) == 0)
                   return(as.data.frame(list()))
                 else
                   return(x)
           })
  if(length(lst) == 0){
    logger.warn("GO LIST is EMPTY")
    return(1)
  }
  #go.xls.list = lapply(lst, fortify, showCategory=Inf)
  go.xls.list <- lapply(names(lst), function(x) lst[[x]]@result)
  names(go.xls.list) = names(go.list)
  WriteXLS(
           go.xls.list,
           file.path(CHARTS_DIR, fxls),
           SheetNames = reformat_sheetnames(names(go.xls.list)))
}

reformat_sheetnames <- function(sheet_names){
  ## []:*?/\
  sheet_names <- stringr::str_replace(sheet_names, "\\[", ".")
  sheet_names <- stringr::str_replace(sheet_names, "\\]", ".")
  sheet_names <- stringr::str_replace(sheet_names, ":", ".")
  sheet_names <- stringr::str_replace(sheet_names, "/", ".")
  sheet_names <- stringr::str_replace(sheet_names, "\\?", ".")
  sheet_names <- stringr::str_replace(sheet_names, "\\\\", ".")
  return(sheet_names)
}


generate_scrna_rawdata <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             data.list = list()
             read_func <- ifelse(CM_FORMAT == "10X", Read10X, Read10X_h5)
             for(i in seq_along(data_src)){
               print(paste(i, Sys.time()))
               key = names(data_src)[i]
               val = data_src[[key]]
               data <- read_func(file.path(data_src[[key]]))
               scrna <- CreateSeuratObject(counts = data, min.cells = MINCELLS, min.features = MINGENES, project = PROJECT)
               name = names(data_src)[i]
               scrna@meta.data[, "name"] <- name
               scrna@meta.data[, "stage"] <- stage_lst[name]
               data.list[[i]] <- scrna
               rm(scrna)
             }
             scrna <- data.list[[1]]
             if(length(data.list) > 1){
               scrna <- merge(x = data.list[[1]], y = unlist(data.list[2:length(data.list)]),
                              merge.data = FALSE, add.cell.ids = names(data_src), project = PROJECT)
             }
             scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-|^MT-")
             scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rpl|^Rps|^RPL|^RPS")
             scrna@tools[["meta_order"]] <- list(name = names(data_src), stage=unique(stage_lst))

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             rm(data.list)
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_ambient_rna <- function(scrna){
  ret_code <- 0

  tryCatch(
           {
            # We convert the input SeuratObject to a SingleCellExperiment object
            # This is the require input for decontX
            assay.used <- DefaultAssay(scrna)
            DefaultAssay(scrna) <- "RNA"
            scrna.sce <- as.SingleCellExperiment(scrna)

            # Now, we estimate and correct the amount of ambient RNA.
            # Since we have not given the function a clustering result,
            # decontX will do the clustering for us.
	    # Since contamination is dependent on the specific experiment, we provide the decontX function with
	    # the batch	information. The contamination is then calculated per sample.
            scrna.decont <- decontX(scrna.sce, batch = scrna$name)

            # We add the estimated contamination and the decontaminated data to the SeuratObject
            scrna[["decontX"]] <- CreateAssayObject(counts = scrna.decont@assays@data$decontXcounts)
            scrna <- AddMetaData(object = scrna, metadata = scrna.decont$decontX_contamination,
                                 col.name = "AmbientRNA")
            scrna <- AddMetaData(object = scrna, metadata = scrna.decont$decontX_clusters,
                                 col.name = "decontX_clusters")
	          DefaultAssay(scrna) <- assay.used

            ## assay to disk
            if (!ALLINONE){
              fname = file.path(SAVE_DIR, "assays", "decontX.Rds")
              # construct list
              assay_info <- list(
                 name = "decontX",
                 assay = scrna[["decontX"]],
                 fname = fname,
                 meta = scrna@meta.data,
                 info = "abundance of ambient RNA")
               #save to disk
              save_object(assay_info, fname, file_format = COMPRESSION_FORMAT)
              if ("assay_info" %ni% names(scrna@tools)){
                 scrna@tools[["assay_info"]] <- list()
              }
              scrna@tools[["assay_info"]][["decontX"]] <- fname
              ## remove from memory
              rm(assay_info)
              scrna[["decontX"]] <- NULL
              gc()
            }

            return(scrna)
           },

            error = function(cond) {
              ret_code <<- -1
              logger.error(cond)
              logger.error(traceback())
           },

            finally={
              return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_phase_singleton <- function(scrna){
  ret_code <- 0
  pconf <- configr::read.config("static/phase.ini")
  pconf <- pconf[["phase_singleton"]]
  pconf <- pconf[pconf== 1]


  for (key in names(pconf)){
    f_name = paste("generate_", key, sep="")
    f_call = paste(f_name, "(scrna)", sep="")
    logger.info(paste("executing", f_name))
    ret_list <-  eval(parse(text=f_call))
    scrna  <- ret_list[[1]]
    ret_code <- ret_list[[2]]
    scrna <- sanity_function(scrna, key)

    if(ret_code != 0){ debug_save_stop(scrna, key) }

    if(key == "scrna_rawdata"){
      scrna@tools$parameter[[cur_time]] <- unlist(pa)
      scrna@tools$execution[[cur_time]] <- conf
      NRds = paste0(key, ".Rds")
      save_object(
        object = scrna,
        file_name = file.path(SAVE_DIR, NRds),
        file_format = COMPRESSION_FORMAT
      )
    }
    logger.info(paste("finished", f_name))
  }
  return(list(scrna, ret_code))
}



generate_scrna_phase_preprocess <- function(scrna){
  ret_code <- 0
  pconf <- configr::read.config("static/phase.ini")
  pconf <- pconf[["phase_preprocess"]]
  pconf <- pconf[pconf== 1]


  for (key in names(pconf)){
    f_name = paste("generate_", key, sep="")
    f_call = paste(f_name, "(scrna)", sep="")
    logger.info(paste("executing", f_name))
    ret_list <-  eval(parse(text=f_call))
    scrna  <- ret_list[[1]]
    ret_code <- ret_list[[2]]
    scrna <- sanity_function(scrna, key)

    if(ret_code != 0){ debug_save_stop(scrna, key) }

    if(key == "scrna_rawdata"){
      scrna@tools$parameter[[cur_time]] <- unlist(pa)
      scrna@tools$execution[[cur_time]] <- conf
      NRds = paste0(key, ".Rds")
      save_object(
        object = scrna,
        file_name = file.path(SAVE_DIR, NRds),
        file_format = COMPRESSION_FORMAT
      )
    }
    gc()
    logger.info(paste("finished", f_name))
  }
  return(list(scrna, ret_code))
}

generate_scrna_phase_clustering <- function(scrna){
  ret_code <- 0
  pconf <- configr::read.config("static/phase.ini")
  pconf <- pconf[["phase_clustering"]]
  pconf <- pconf[pconf== 1]

  for (key in names(pconf)){
      f_name = paste("generate_", key, sep="")
      f_call = paste(f_name, "(scrna)", sep="")
      logger.info(paste("executing", f_name))
      ret_list <-  eval(parse(text=f_call))
      scrna  <- ret_list[[1]]
      ret_code <- ret_list[[2]]
      scrna <- sanity_function(scrna, key)

      if(ret_code != 0){ debug_save_stop(scrna, key) }

      gc()
      logger.info(paste("finished", f_name))
   }
  return(list(scrna, ret_code))
}


generate_scrna_phase_existed_clusters <- function(scrna){
  ret_code <- 0
  pconf <- configr::read.config("static/phase.ini")
  pconf <- pconf[["phase_existed_clusters"]]
  pconf <- pconf[pconf== 1]

  if(exists("existed_cluster_slots_map")){
      nms <- names(existed_cluster_slots_map)
      if (length(nms) > 0){
          metas <- nms[startsWith(nms, "meta.")]
          reductions <- nms[startsWith(nms, "reduction.")]
          for(meta in metas){
              assertthat::assert_that(existed_cluster_slots_map[meta] %in% names(scrna@meta.data))
              meta_slot <- stringr::str_sub(meta, start=(nchar("meta.")+1))
              scrna@meta.data[, meta_slot] <- scrna@meta.data[, existed_cluster_slots_map[meta]]
          }
          for(reduction in reductions){
              assertthat::assert_that(existed_cluster_slots_map[reduction] %in% names(scrna@reductions))
              reduction_slot <- stringr::str_sub(reduction, start=(nchar("reduction.")+1))
              scrna[[reduction_slot]] <- scrna[[existed_cluster_slots_map[reduction]]]
          }
      }

  }

  for (key in names(pconf)){
      f_name = paste("generate_", key, sep="")
      f_call = paste(f_name, "(scrna)", sep="")
      logger.info(paste("executing", f_name))
      ret_list <-  eval(parse(text=f_call))
      scrna  <- ret_list[[1]]
      ret_code <- ret_list[[2]]
      scrna <- sanity_function(scrna, key)

      if(ret_code != 0){ debug_save_stop(scrna, key) }
      gc()
      logger.info(paste("finished", f_name))
   }
  return(list(scrna, ret_code))
}


generate_scrna_phase_comparing <- function(scrna){
  ret_code <- 0
  pconf <- configr::read.config("static/phase.ini")
  pconf <- pconf[["phase_comparing"]]
  pconf <- pconf[pconf== 1]
  for (key in names(pconf)){
      f_name = paste("generate_", key, sep="")
      f_call = paste(f_name, "(scrna)", sep="")
      logger.info(paste("executing", f_name))
      ret_list <-  eval(parse(text=f_call))
      scrna  <- ret_list[[1]]
      ret_code <- ret_list[[2]]
      scrna <- sanity_function(scrna, key)

      if(ret_code != 0){ debug_save_stop(scrna, key) }
      gc()
      logger.info(paste("finished", f_name))
   }
  return(list(scrna, ret_code))

}



generate_scrna_filter <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             scrna <- subset(scrna, subset = nFeature_RNA > nFeatureRNAfloor &
                             nFeature_RNA < nFeatureRNAceiling&
                             nCount_RNA > nCountRNAfloor&
                             nCount_RNA < nCountRNAceiling &
                             percent.mt > pct_mitofloor&
                             percent.mt < pct_mitoceiling &
                             percent.ribo > pct_ribofloor &
                             percent.ribo < pct_riboceiling)
           },
            error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

## !!!!Dangerous, in this case, once deleted, never recovered
generate_scrna_del_mitogenes <- function(scrna){
  return_code = 0
  MT_PATTERN = "^mt-"
  if(SPECIES == "Human"){ ## Default is Mouse
    MT_PATTERN = "^MT-"
  }
  allgenes <- rownames(scrna)
  mitogenes <-   allgenes[grepl(MT_PATTERN, allgenes)]
  restgenes <- setdiff(allgenes, mitogenes)
  scrna <- scrna[restgenes, ]

  return(list( scrna, return_code))
}

generate_scrna_preprocess <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-|^MT-")
             scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rpl|^Rps|^RPL|^RPS")

             mt.genes <- grep(pattern = "^mt-|^MT-", x = rownames(x = scrna), value = TRUE)
             ribo.genes <- grep("^Rpl|^Rps|^RPL|^RPS",  x = rownames(x = scrna), value = TRUE)

             scrna[["percent.exclude"]]  <- PercentageFeatureSet(scrna, features = c(mt.genes, ribo.genes))

             scrna <- NormalizeData(scrna)
             scrna <- FindVariableFeatures(scrna, selection.method = "vst")
             scrna <- ScaleData(scrna, features = rownames(scrna))


             all.genes <- rownames(x = scrna)
             s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
             s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
             g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
             g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
             scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes), reduction.name="RAW_PCA")
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

#pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE) ## another way to preprocess the data

generate_scrna_cellcycle <-function(scrna){
  ret_code = 0
  tryCatch(
           {

             all.genes <- rownames(x = scrna)

             s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
             s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]

             g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
             g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]

             scrna <- NormalizeData(scrna)
             scrna <- ScaleData(scrna, features = rownames(scrna))
             scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes), reduction.name="BCELLCYCLE_PCA")
             scrna <- CellCycleScoring(scrna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
             scrna$G1.Score = 1 - scrna$S.Score - scrna$G2M.Score
             scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
             if("percent.mt" %ni% names(scrna@meta.data)){ ## for existed clusters analysis
                scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-|^MT-")
             }
             if("percent.ribo" %ni% names(scrna@meta.data)){
                scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rpl|^Rps|^RPL|^RPS")
             }

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_cycleRegressOut <- function(scrna){
  ret_code = 0
  tryCatch(
           {

             all.genes <- rownames(x = scrna)

             s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
             s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]

             g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
             g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]


             ### Scale data add exclude
             scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
             scrna <- ScaleData(scrna, vars.to.regress = c("G2M.Score", "S.Score"), features = rownames(scrna))
             scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes), nfeatures.print = 10, reduction.name="CELLCYCLED_PCA")
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })

  return(list(scrna, ret_code))
}

## Regressout cell cycle and mitochondrial
generate_scrna_regressOut<- function(scrna){
  ret_code = 0
  tryCatch(
           {
             scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA", get_regressout_vector()), features = rownames(scrna))
             scrna <- RunPCA(scrna, features = VariableFeatures(scrna), nfeatures.print = 10, reduction.name="RegressOut_PCA")
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}




select_features <- function(data.list){
  ret_code = 0
  tryCatch(
           {
             data.list <- SplitObject(scrna, split.by = "name")
             data.list <- lapply(X = data.list, FUN = function(x) {
                                   x <- NormalizeData(x, verbose = FALSE)
                                   x <- FindVariableFeatures(x, verbose = FALSE)})
             features <- SelectIntegrationFeatures(object.list = data.list)

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(features, ret_code))
           })
  return(features)
}


scale_pca <- function(features){
  ret_code = 0
  tryCatch(
           {
             data.list <- lapply(X = data.list, FUN = function(x) {
                                   x <- ScaleData(x, features = features, verbose = FALSE)
                                   x <- RunPCA(x, features = features, verbose = FALSE)
                                   })
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(data.list)
           })
  return(data.list)
}


generate_scrna_integration_harmony <- function(scrna){
  ret_code = 0
  tryCatch(
           {
            largestDim=ncol(Seurat::Embeddings(scrna[["RegressOut_PCA"]]))
            keep_harmony_dims = HARMONY_DIM[HARMONY_DIM<=largestDim]
            scrna <- harmony::RunHarmony(scrna,  "name", plot_convergence = "True", reduction = "RegressOut_PCA")
            scrna <- RunUMAP(scrna, reduction = "harmony", dims = keep_harmony_dims, reduction.name= "harmony_UMAP")
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })

  return(list(scrna, ret_code))
}

generate_scrna_integration_seurat <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             data.list <- SplitObject(scrna, split.by = "name")
             data.list <- lapply(X = data.list, FUN = function(x) {
                                   x <- NormalizeData(x)
                                   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

             ## scale = False to use previous scaled data
	     ## If the number of cells is < 200 for a sample, we need to reduce the default values.
	     ## If the number of cells is < 200 for any sample, we set the k.filter and k.weight values to this new minimum.
	     ## This way, we can still integrate very small data sets.
	     k.filter <- min(table(scrna$name))
             k.filter <- ifelse(k.filter < 200, k.filter, 200)
             anchors <- FindIntegrationAnchors(object.list = data.list, dims = INTEGRATED_DIM, scale=F,
					       k.filter = k.filter)## THIS IS CCA DIMENSIONS
             scrna_inte <- IntegrateData(anchorset = anchors, dims = INTEGRATED_DIM, k.weight = k.filter) ## THIS IS CCA DIMENSION
             ## keep the order of integration obj
             scrna_inte <- scrna_inte[, colnames(scrna)]
             scrna[['integrated']] <- scrna_inte[['integrated']]
             scrna@commands <- c(scrna@commands, scrna_inte@commands)
             scrna@tools <- c(scrna@tools, scrna_inte@tools)
             DefaultAssay(scrna) <- "integrated"
             scrna <- ScaleData(scrna, verbose = FALSE)
             scrna <- RunPCA(scrna, npcs = max(50, max(FINDNEIGHBORS_DIM)), verbose = FALSE, reduction.name="INTE_PCA")
             scrna@reductions$INTE_PCA@assay.used <- "RNA"
             DefaultAssay(scrna) <- "RNA"
             scrna <- RunUMAP(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM, reduction.name="INTE_UMAP")
             rm(scrna_inte)
             rm(data.list)
             ## assay to disk
             if (!ALLINONE){
               fname = file.path(SAVE_DIR, "assays", "integrated.Rds")
               # construct list
               assay_info <- list(
                 name = "integrated",
                 assay = scrna[["integrated"]],
                 fname = fname,
                 meta = scrna@meta.data,
                 info = "seurat integrated assay")
               #save to disk
               save_object(assay_info, fname, file_format = COMPRESSION_FORMAT)
               if ("assay_info" %ni% names(scrna@tools)){
                  scrna@tools[["assay_info"]] <- list()
               }
               scrna@tools[["assay_info"]][["integrated"]] <- fname
               ## remove from memory
               rm(assay_info)
               scrna[["integrated"]] <- NULL
               gc()
             }


           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_ScaleSingleton <- function(scrna){
  ret_code = 0
  tryCatch(
           {

             # Run the standard workflow for visualization and clustering
             scrna <- ScaleData(scrna, verbose = FALSE)
             scrna <- RunPCA(scrna, npcs = max(50, max(FINDNEIGHBORS_DIM)), verbose = FALSE, reduction.name="SINGLE_PCA")
             scrna <- RunUMAP(scrna, reduction = "SINGLE_PCA", dims = FINDNEIGHBORS_DIM, reduction.name="SINGLE_UMAP")
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}


generate_scrna_sltn_batch_clustering <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             scrna <- FindNeighbors(scrna, reduction = "SINGLE_PCA", dims = FINDNEIGHBORS_DIM)
             scrna <- FindClusters(scrna, resolution = CLUSTER_RESOLUTION_RANGE) ##

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}



generate_scrna_singleton_clustering <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             a_meta <- sprintf("RNA_snn_res.%s", CLUSTER_RESOLUTION)
             if(a_meta %in% colnames(scrna@meta.data)){
                scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- scrna@meta.data[, a_meta]
             }else{
               scrna <- FindNeighbors(scrna, reduction = "SINGLE_PCA", dims = FINDNEIGHBORS_DIM)
               scrna <- FindClusters(scrna, resolution = CLUSTER_RESOLUTION) ##
             }
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}
generate_scrna_clustering <- function(scrna){
  ret_code = 0
  tryCatch(
           {
               #DefaultAssay(scrna) <- "integrated"
               DefaultAssay(scrna) <- "RNA"
               scrna <- FindNeighbors(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM, graph.name='integrated_snn') %>%
                            FindClusters(resolution = CLUSTER_RESOLUTION, graph.name='integrated_snn') ##

               scrna$seurat_inte_clusters <- scrna$seurat_clusters

               DefaultAssay(scrna) <- "RNA"
               scrna = FindNeighbors(scrna, reduction = "harmony", dims = HARMONY_DIM) %>%
                            FindClusters(resolution = CLUSTER_RESOLUTION) ##
               scrna$harmony_inte_clusters <- scrna$seurat_clusters

             if(INTEGRATION_OPTION == "harmony"){
                scrna$seurat_clusters <- scrna$harmony_inte_clusters
                scrna[["DEFAULT_UMAP"]] <- scrna[["harmony_UMAP"]]
             }else if(INTEGRATION_OPTION == "seurat"){
                scrna$seurat_clusters <- scrna$seurat_inte_clusters
                scrna[["DEFAULT_UMAP"]] <- scrna[["INTE_UMAP"]]
             }

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_clusterwise_xcell<- function(scrna){
  #remove cells of each cluster according distinct criterion
  ret_code = 0
  Idents(scrna) <- DEFUALT_CLUSTER_NAME
  tryCatch(
           {
             df <- do.call(rbind, scrna_clusterwise_filtercell_settings)
             ### clusters should be included in all clusters
             assertthat::assert_that(all(df$cluster %in% unique(scrna@meta.data[, DEFUALT_CLUSTER_NAME])))

             del_cells <- c()
             for(x in rownames(df)){
               ## get the cells that in a cluster && meet metrics too
               meta <- scrna@meta.data
               cluster_bools <- meta[, DEFUALT_CLUSTER_NAME] == df[x, ]$cluster
               meta <- meta[cluster_bools, ]

               ##DEFAULT MITO
               min_pct_bools <- meta$percent.mt >= df[x, ]$min_pct
               max_pct_bools <- meta$percent.mt <= df[x, ]$max_pct
               if(df[x, ]$type == "ribo"){
                 min_pct_bools <- meta$percent.ribo >= df[x, ]$min_pct
                 max_pct_bools <- meta$percent.ribo <= df[x, ]$max_pct
               }
               submeta <- meta[min_pct_bools & max_pct_bools, ]
               sub_dell_cells <- setdiff(rownames(meta), rownames(submeta))

               del_cells <- c(del_cells, sub_dell_cells)
             }

             keep_cells <-  setdiff(colnames(scrna), del_cells)
             scrna <- subset(scrna, cells = keep_cells)
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}



generate_scrna_batchclustering <- function(scrna){
  ret_code = 0
  tryCatch(
           {

             #DefaultAssay(scrna) <- "integrated"
             DefaultAssay(scrna) <- "RNA"
             scrna <- FindNeighbors(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM, graph.name='integrated_snn') %>%
                                  FindClusters(resolution = CLUSTER_RESOLUTION_RANGE, graph.name='integrated_snn') ##

             DefaultAssay(scrna) <- "RNA"
             scrna <-  FindNeighbors(scrna, reduction = "harmony", dims = HARMONY_DIM) %>%
                                  FindClusters(resolution = CLUSTER_RESOLUTION_RANGE) ##

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}




generate_scrna_fishertest_inte_clusters <- function(scrna){
  ret_code = 0
  CLUSTER_TO_TEST_VEC <- c("seurat_inte_clusters",  "harmony_inte_clusters")
  for(CLUSTER_TO_TEST in CLUSTER_TO_TEST_VEC){

     STAGE_TO_TEST <- "stage"
     stages <- scrna@tools[["meta_order"]][[STAGE_TO_TEST]]

     ## to avoid matrix still a table
     count.matrix <- as.matrix(Matrix(table(scrna@meta.data[, CLUSTER_TO_TEST], scrna@meta.data[, STAGE_TO_TEST])))
     count.matrix <- count.matrix[rowSums(count.matrix)>0, ]
     df <- as.data.frame(count.matrix, stringsAsFactors=F)
     help_sort_func <- ifelse(all.is.numeric(names(df)), as.numeric, function(x){x})
     inm <-help_sort_func(rownames(df))
     if(all.is.numeric(names(df))){
       df$Cluster <- factor(rownames(df), levels=(min(inm)):(max(inm)))
     }else{
       df$Cluster <- factor(rownames(df))#, levels=(min(inm)):(max(inm)))
     }


     stages_comb <- list_stages()
     for (x in 1:length(stages_comb)){
       a <- stages_comb[[x]][1]
       b <- stages_comb[[x]][2]
       df[, glue("o.{a}")] <- sum(count.matrix[, a]) - df[, a]
       df[, glue("o.{b}")] <- sum(count.matrix[, b]) - df[, b]
       for(cluster in df$Cluster){
         mtx <- matrix(unlist(df[cluster, c(a, b, glue("o.{a}"), glue("o.{b}"))]), nrow=2)
         ft <- fisher.test(mtx,workspace=1e9)
         df[cluster, glue("pval_{a}.vs.{b}")] <- ft$p.value
         df[cluster, glue("odds.ratio_{a}.vs.{b}")] <- ft$estimate
       }
       df[, glue("pval.adjust_{a}.vs.{b}")] <- p.adjust(df[, glue("pval_{a}.vs.{b}")],  method = "bonferroni")
     }

    if(!ALLINONE){
      fname = file.path(SAVE_DIR, "partition", glue("fishertest_{CLUSTER_TO_TEST}.Rds"))
       save_object(df,
                   file_name = fname,
                   file_format = COMPRESSION_FORMAT)
       scrna@tools[[sprintf("fishertest_%s", CLUSTER_TO_TEST)]] <- fname
    }else{
       scrna@tools[[sprintf("fishertest_%s", CLUSTER_TO_TEST)]] <- df
    }
  }
  rm(count.matrix)
  return(list(scrna, ret_code))
}


generate_scrna_fishertest_clusters <- function(scrna){
  ret_code = 0
  CLUSTER_TO_TEST <- DEFUALT_CLUSTER_NAME
  STAGE_TO_TEST <- "stage"
  stages <- scrna@tools[["meta_order"]][[STAGE_TO_TEST]]
  assertthat::assert_that(CLUSTER_TO_TEST %in% names(scrna@meta.data))
  assertthat::assert_that(STAGE_TO_TEST %in% names(scrna@meta.data))


  ## to avoid matrix still a table
  count.matrix <- as.matrix(Matrix(table(scrna@meta.data[, CLUSTER_TO_TEST], scrna@meta.data[, STAGE_TO_TEST])))
  count.matrix <- count.matrix[rowSums(count.matrix)>0, ]
  df <- as.data.frame(count.matrix, stringsAsFactors=F)
  help_sort_func <- ifelse(all.is.numeric(names(df)), as.numeric, function(x){x})
  inm <-help_sort_func(rownames(df))
  if(all.is.numeric(names(df))){
    df$Cluster <- factor(rownames(df), levels=(min(inm)):(max(inm)))
  }else{
    df$Cluster <- factor(rownames(df))#, levels=(min(inm)):(max(inm)))
  }


  stages_comb <- list_stages()
  for (x in 1:length(stages_comb)){
    a <- stages_comb[[x]][1]
    b <- stages_comb[[x]][2]
    df[, glue("o.{a}")] <- sum(count.matrix[, a]) - df[, a]
    df[, glue("o.{b}")] <- sum(count.matrix[, b]) - df[, b]
    for(cluster in df$Cluster){
      mtx <- matrix(unlist(df[cluster, c(a, b, glue("o.{a}"), glue("o.{b}"))]), nrow=2)
      ft <- fisher.test(mtx,workspace=1e9)
      df[cluster, glue("pval_{a}.vs.{b}")] <- ft$p.value
      df[cluster, glue("odds.ratio_{a}.vs.{b}")] <- ft$estimate
    }
    df[, glue("pval.adjust_{a}.vs.{b}")] <- p.adjust(df[, glue("pval_{a}.vs.{b}")],  method = "bonferroni")
  }

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", glue("fishertest_{CLUSTER_TO_TEST}.Rds"))
     save_object(df,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("fishertest_%s", CLUSTER_TO_TEST)]] <- fname
  }else{
     scrna@tools[[sprintf("fishertest_%s", CLUSTER_TO_TEST)]] <- df
  }
  rm(count.matrix)
  return(list(scrna, ret_code))

}
generate_scrna_proptest_clusters <- function(scrna){
  ret_code = 0
  source("R/scProportion.R")

  CLUSTER_TO_TEST <- DEFUALT_CLUSTER_NAME
  stages_comb <- list_stages()
  prop_test_all <- sc_utils(scrna)
  proptest_list <- list()
  for (x in stages_comb){
    prop_test <- prop_test_all
    prop_test@meta_data <- prop_test@meta_data %>% dplyr::filter(stage %in% x)
    prop_test <- permutation_test(
       prop_test, cluster_identity = CLUSTER_TO_TEST,
       sample_1 = x[1], sample_2 = x[2],
       sample_identity = "stage")
    proptest_list[[glue("{x[1]}.vs.{x[2]}")]] <- data.table::as.data.table(prop_test@results$permutation)
  }

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", glue("proptest_{DEFUALT_CLUSTER_NAME}.Rds"))
     save_object(proptest_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("proptest_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("proptest_%s", DEFUALT_CLUSTER_NAME)]] <- proptest_list
  }
  rm(proptest_list)
  rm(prop_test_all)
  return(list(scrna, ret_code))
}

generate_scrna_merge_clusters <- function(scrna){
  ret_code = 0
  scrna$merged_clusters <- as.character(scrna$seurat_clusters)
  for (nm in names(scrna_merge_clusters)){
    scrna$merged_clusters[scrna$merged_clusters %in% as.character(scrna_merge_clusters[[nm]])] <- nm
  }
  return(list(scrna, ret_code))
}




generate_scrna_remove_clusters <- function(scrna){
  ret_code = 0
  Idents(scrna) <- "seurat_clusters"
  keeps <- setdiff(unique(scrna$seurat_clusters), scrna_remove_clusters)
  scrna <- subset(scrna, idents= keeps)
  scrna$removed_clusters <- scrna$seurat_clusters
  if(is.factor(scrna$removed_clusters)){
    scrna$removed_clusters = droplevels(scrna$removed_clusters)
  }
  return(list(scrna, ret_code))
}


generate_scrna_cluster_annotation <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             clusters <- unique(as.character(scrna@meta.data[, from_cluster_slot]))
             nms <- names(cluster_annotation)
             if(length(setdiff(clusters, nms)) > 0){
               ret_code <<- 1
               stop("Clusters are not consistent with annotation clusters")
             }

             annotation <- plyr::mapvalues(scrna@meta.data[, from_cluster_slot],
                                           names(cluster_annotation),
                                           cluster_annotation)
             scrna$annotation <- annotation
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_remove_recluster <- function(scrna){
  ret_code = 0
  Idents(scrna) <- "seurat_clusters"
  keeps <- setdiff(unique(scrna$seurat_clusters), scrna_remove_clusters)
  scrna <- subset(scrna, idents= keeps)

  execution_vec = c(
        "generate_scrna_integration_seurat(scrna)",
        "generate_scrna_integration_harmony(scrna)",
        "generate_scrna_clustering(scrna)"
  )
  for(func_call in execution_vec){
    logger.info(paste("executing ", func_call))
    ret_list <- eval(parse(text=func_call))
    scrna <- ret_list[[1]]
    ret_code <- ret_list[[2]]
    if(ret_code != 0){
      return(list(scrna, ret_code))
    }
  }
  scrna$remove_recluster <- scrna$seurat_clusters
  if(is.factor(scrna$remove_recluster)){
    scrna$remove_recluster = droplevels(scrna$remove_recluster)
  }
  ## recluster, need to redo everything related to clusterwise
  execution_vec = c(
        "generate_scrna_fishertest_clusters(scrna)",
        "generate_scrna_MCAannotate(scrna)",
        "generate_scrna_ExternalAnnotation(scrna)",
        "generate_scrna_HCLannotate(scrna)"
  )
  for(func_call in execution_vec){
    logger.info(paste("executing ", func_call))
    ret_list <- eval(parse(text=func_call))
    scrna <- ret_list[[1]]
    ret_code <- ret_list[[2]]
    if(ret_code != 0){
      return(list(scrna, ret_code))
    }
  }
  return(list(scrna, ret_code))
}


generate_scrna_MCAannotate <- function(scrna){
  ret_code = 0
  suppressPackageStartupMessages(require(scMCA))
  mca_result <- scMCA(GetAssayData(object=scrna, slot="counts"), numbers_plot = 3)
  #pattern = paste0("(", gsub("-", "_", MCA_NAME), ")")
  pattern = paste0("\\(", MCA_NAME, "\\)")
  corr=mca_result$cors_matrix[grep(pattern,rownames(mca_result$cors_matrix)),]
  if(dim(corr)[1] == 0){
    logger.error("Cannot find MCA name, please check!")
    stop("exit 1")
  }
  res=rownames(corr)[max.col(t(corr))]
  scrna$MCA_annotate <- res
  Idents(object=scrna) <- "name"
  rm(mca_result)
  return(list(scrna, ret_code))
}

generate_scrna_HCLannotate <- function(scrna){
  ret_code = 0
  suppressPackageStartupMessages(require(scHCL))
  hcl_result <- scHCL(GetAssayData(object=scrna, slot="counts"), numbers_plot = 3)
  pattern = gsub("-",".",HCL_NAME)
  corr=hcl_result$cors_matrix[grep(pattern,rownames(hcl_result$cors_matrix),fixed=TRUE),]
  rnms <- gsub("_[a-zA-Z]+\\.$", "", rownames(corr))
  rnms <- gsub(paste0("\\.", pattern, "\\."), "", rnms)
  rnms <- gsub("_.*high", "", rnms)
  rnms <- gsub("[0-9]\\.$", "", rnms)
  rnms <- gsub("\\.[0-9]$", "", rnms)
  rnms <- gsub("[0-9]\\.$", "", rnms)
  rnms <- gsub("\\.\\.", ".", rnms)
  rnms <- gsub("\\.$", "", rnms)
  rownames(corr) <- rnms
  if(dim(corr)[1] == 0){
    logger.error("Cannot find HCL name, please check!")
    stop("exit 1")
  }
  res=rownames(corr)[max.col(t(corr))]
  scrna$HCL_annotate <- res
  Idents(object=scrna) <- "name"
  return(list(scrna, ret_code))
}

generate_scrna_MAGIC <- function(scrna){
  ret_code = 0
  DefaultAssay(scrna) <- "RNA"
  all_genes <- rownames(scrna)
  scrna <- magic(scrna, genes=all_genes)
  rm(all_genes)
  ## assay to disk
  if (!ALLINONE){
    fname = file.path(SAVE_DIR, "assays", "MAGIC_RNA.Rds")
    # construct list
    assay_info <- list(
      name = "MAGIC_RNA",
      assay = scrna[["MAGIC_RNA"]],
      fname = fname,
      meta = scrna@meta.data,
      info = "MAGIC imputed RNA assay")
    #save to disk
    save_object(assay_info, fname, file_format = COMPRESSION_FORMAT)
    if ("assay_info" %ni% names(scrna@tools)){
       scrna@tools[["assay_info"]] <- list()
    }
    scrna@tools[["assay_info"]][["MAGIC_RNA"]] <- fname
    ## remove from memory
    rm(assay_info)
    scrna[["MAGIC_RNA"]] <- NULL
    gc()
  }
  return(list(scrna, ret_code))
}


generate_scrna_ExternalAnnotation <- function(scrna){
  ret_code = 0
  df <- read.csv(file=ANNOTATION_EXTERNAL_FILE, sep="\t", stringsAsFactors=F)
  df <- df[!apply(df, 1, function(x) all(x=="")), ]

  if(!(ORGAN %in% df$Tissue.of.Origin)){
    logger.error("Cannot find the organ data from the external annotation data")
    stop("exit 1")
  }

  if(!(SPECIES %in% c("Human", "Mouse"))){
    logger.error("Cannot find the species data from the external annotation data")
    stop("exit 1")
  }

  dfs <- split(df, df$Tissue.of.Origin)

  DefaultAssay(scrna) <- "RNA"
  mtx <- GetAssayData(object = scrna, slot = "data")
  anno_genes <- unique(dfs[[ORGAN]][, sprintf("%s.Gene", SPECIES)])
  use_genes <- intersect(anno_genes, rownames(mtx))
  df <- dfs[[ORGAN]]

  celltype.list <- foreach(a_cell = colnames(mtx)) %dopar% {
    score <- mtx[use_genes, a_cell]
    score <- score[score>0]
    if(length(score) == 0){
      return("Unknown")
    }
    sdf <- df[df[,glue("{SPECIES}.Gene")] %in% names(score), c(glue("{SPECIES}.Gene"), "Cell.Type")]
    sdf$score <- score[sdf[,glue("{SPECIES}.Gene")]]
    itype <- aggregate(sdf$score, by=list(CellType=sdf$Cell.Type),
                       FUN=function(x){return(mean(x)*log(length(x)+1))})
    celltype <- itype[which.max(itype$x), ]$CellType
    return(celltype)
  }
  celltypes <- unlist(celltype.list)
  names(celltypes)<-colnames(mtx)
  scrna$external_annotation <- unlist(celltypes)
  rm(mtx)
  return(list(scrna, ret_code))

}


generate_scrna_markergenes <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             DefaultAssay(scrna) <- "RNA"
             Idents(scrna) <- DEFUALT_CLUSTER_NAME
             de.df = FindAllMarkers(scrna, logfc.threshold=0)

             cluster.de <- split(de.df, de.df$cluster)


             if(!ALLINONE){
               fname = file.path(SAVE_DIR, "partition", sprintf("de_%s.Rds", DEFUALT_CLUSTER_NAME))
               save_object(cluster.de,
                            file_name = fname,
                            file_format = COMPRESSION_FORMAT)
                scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]] <- fname
             }else{
                scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]] <- cluster.de
             }
             flist <- lapply(cluster.de, subset, subset = p_val_adj < 0.05)
             flist <-  flist[sapply(flist, function(m) nrow(m) >0)]

             WriteXLS(
                      flist,
                      file.path(CHARTS_DIR, sprintf("de_%s.xlsx", DEFUALT_CLUSTER_NAME)),
                      SheetNames = reformat_sheetnames(names(flist)))
             rm(cluster.de)
             rm(flist)
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_batch_markergenes <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             len <- length(CLUSTER_RESOLUTION_RANGE)
             cluster.de.list <- vector("list", length = len)
             names(cluster.de.list) <- as.character(CLUSTER_RESOLUTION_RANGE)
             DefaultAssay(scrna) <- "RNA"
             cluster.de.list <- foreach(i=CLUSTER_RESOLUTION_RANGE) %do%{
               cluster_name <- sprintf("integrated_snn_res.%s", i)
               if(INTEGRATION_OPTION == "harmony"){
                 cluster_name <- sprintf("RNA_snn_res.%s", i)
               }
               Idents(scrna) <- cluster_name
               de.df = FindAllMarkers(scrna, logfc.threshold=0)
               cluster.de <- split(de.df, de.df$cluster)
             }

             if(!ALLINONE){
               fname = file.path(SAVE_DIR, "partition", "de_batch.Rds")
               save_object(cluster.de.list,
                           file_name = fname,
                           file_format = COMPRESSION_FORMAT)
               scrna@tools[["de_batch"]] <- fname
             }else{
               scrna@tools[["de_batch"]] <- cluster.de.list
             }
             rm(cluster.de.list)
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_singleton_markergenes <- function(scrna){
  ret_code = 0
  tryCatch(
           {
             len <- length(CLUSTER_RESOLUTION_RANGE)
             cluster.de.list <- vector("list", length = len)
             names(cluster.de.list) <- as.character(CLUSTER_RESOLUTION_RANGE)
             cluster.de.list <- foreach(i=CLUSTER_RESOLUTION_RANGE) %do%{
               DefaultAssay(scrna) <- "RNA"
               cluster_name <- sprintf("RNA_snn_res.%s", i)
               Idents(scrna) <- cluster_name
               de.df = FindAllMarkers(scrna, logfc.threshold=0)
               cluster.de <- split(de.df, de.df$cluster)
             }


             if(!ALLINONE){
               fname = file.path(SAVE_DIR, "partition", "de_batch.Rds")
               save_object(cluster.de.list,
                            file_name = fname,
                            file_format = COMPRESSION_FORMAT)
                scrna@tools[["de_batch"]] <- fname
             }else{
                scrna@tools[["de_batch"]] <- cluster.de.list
             }
             rm(cluster.de)

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}



generate_scrna_genesorteR <- function(scrna){
  suppressPackageStartupMessages(library(genesorteR))
  ret_code = 0
  tryCatch(
           {
             DefaultAssay(scrna) <- "RNA"
             Idents(scrna) <- DEFUALT_CLUSTER_NAME
             genesorter = sortGenes(GetAssayData(scrna, slot="counts"), Idents(scrna))

             pv <- getPValues(genesorter)
             tbl = getTable(genesorter, pv, fc_cutoff = 0, adjpval_cutoff = 0.05,
                            islog = TRUE, pseudocount = 1)
             tbl$gene <- rownames(tbl)
             flist <- split(tbl, tbl$Cluster)
             flist <-  flist[sapply(flist, function(m) nrow(m) >0)]

             WriteXLS(
                      flist,
                      file.path(CHARTS_DIR, sprintf("genesorteR_%s.xlsx", DEFUALT_CLUSTER_NAME)),
                      SheetNames = reformat_sheetnames(names(flist)))

             condGeneProb <- as.data.frame(genesorter$condGeneProb)
             condGeneProb <- cbind(rownames(condGeneProb), condGeneProb)
             colnames(condGeneProb)[1] <- "gene"

             WriteXLS(condGeneProb,
                      file.path(CHARTS_DIR, glue("genesorteR_condGeneProb_{DEFUALT_CLUSTER_NAME}.xlsx")))

             store_list <- list(genesorter, flist, condGeneProb)
             names(store_list) <- c("obj", "tbl", "condprob")


             if(!ALLINONE){
               fname = file.path(SAVE_DIR, "partition", sprintf("genesorteR_%s.Rds", DEFUALT_CLUSTER_NAME))
               save_object(store_list,
                            file_name = fname,
                            file_format = COMPRESSION_FORMAT)
                scrna@tools[[sprintf("genesorteR_%s", DEFUALT_CLUSTER_NAME)]] <- fname
             }else{
                scrna@tools[[sprintf("genesorteR_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
             }
             rm(store_list)

           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },

           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}



generate_scrna_dego_name <- function(scrna){
  ret_code = 0
  lst <- list_1vs1()
  all_de_list = list()
  all_goup_list = list()
  all_godown_list = list()
  all_de_list <- mclapply(lst, function(ident.use){
                            Idents(scrna) <- "name"
                            logger.info("****processing %s vs %s", ident.use[1], ident.use[2])
                            used_idents <- which(scrna$name %in% ident.use)
                            a_sub <- subset(scrna, cells= used_idents)
                            Idents(a_sub) <- DEFUALT_CLUSTER_NAME
                            de.list <- list()
                            de.list = mclapply(unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME])), function(clusters){
                                                 logger.info("processing cluster %s", as.character(clusters))
                                                 sub_idents = which(a_sub@meta.data[, DEFUALT_CLUSTER_NAME] == clusters)
                                                 if(length(sub_idents) == 0){
                                                   return(list(NULL))
                                                   next
                                                 }
                                                 a_sub.subset <- subset(a_sub, cells = sub_idents)
                                                 Idents(a_sub.subset) <- "name"

                                                 cluster.inside.de <- list(NULL)
                                                 tryCatch({
                                                   cluster.inside.de <- FindMarkers(a_sub.subset, ident.1= ident.use[1], logfc.threshold=0)
                                                 }, error=function(cond) {
                                                   logger.warn(cond)
                                                   cluster.inside.de <- list(NULL)
                                                 },finally={
                                                   if (length(cluster.inside.de) != 0){
                                                     cluster.inside.de$gene <- rownames(cluster.inside.de)
                                                     return(cluster.inside.de)
                                                   }else{
                                                     return(list(NULL))
                                                   }
                                                 })


                      }, mc.cores=1)
                            names(de.list) <- unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME]))
                            de.list <- de.list[sapply(de.list, function(x){(!is.null(x[[1]]))})]
                            nm <- paste0(ident.use[1], ".vs.", ident.use[2])
                            return(de.list)
           }, mc.cores=1)

  names(all_de_list) <- sapply(lst, function(x) paste0(x[1], ".vs.", x[2]) )

  for(nm in names(all_de_list)) {
    logger.info("****processing go up&down %s", nm)
    de.list <- all_de_list[[nm]]
    go_ups <- get_go_up(de.list)
    all_goup_list[[nm]] <- go_ups
    go_downs <- get_go_down(de.list)
    all_godown_list[[nm]] <- go_downs
    dego_dump(nm, de.list, go_ups, go_downs)
    rm(de.list)
    rm(go_ups)
    rm(go_downs)
  }

  store_list <- list(all_de_list, all_goup_list, all_godown_list)
  names(store_list) <- c("de", "goup", "godown")

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("dego_name_%s.Rds", DEFUALT_CLUSTER_NAME))
    save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("dego_name_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("dego_name_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }
  rm(store_list)
  rm(all_de_list)
  rm(all_goup_list)
  rm(all_godown_list)


  return(list(scrna, ret_code))
}


generate_scrna_dego_stage <- function(scrna){
  ret_code = 0
  lst <- list_stages()
  all_de_list = list()
  all_goup_list = list()
  all_godown_list = list()
  all_de_list <- mclapply(lst, function(ident.use){
                            Idents(scrna) <- "stage"
                            logger.info("****processing %s vs %s", ident.use[1], ident.use[2])
                            used_idents <- which(scrna$stage %in% ident.use)
                            a_sub <- subset(scrna, cells= used_idents)
                            Idents(a_sub) <- DEFUALT_CLUSTER_NAME
                            de.list <- list()
                            de.list = mclapply(unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME])), function(clusters){
                                                 logger.info("processing cluster %s", as.character(clusters))
                                                 sub_idents = which(a_sub@meta.data[, DEFUALT_CLUSTER_NAME] == clusters)
                                                 if(length(sub_idents) == 0){
                                                   return(list(NULL))
                                                   next
                                                 }
                                                 a_sub.subset <- subset(a_sub, cells = sub_idents)
                                                 Idents(a_sub.subset) <- "stage"

                                                 cluster.inside.de <- list(NULL)
                                                 tryCatch({
                                                   cluster.inside.de <- FindMarkers(a_sub.subset, ident.1= ident.use[1], logfc.threshold=0)
                                                 }, error=function(cond) {
                                                   logger.warn(cond)
                                                   cluster.inside.de <- list(NULL)
                                                 },finally={
                                                   if (length(cluster.inside.de) != 0){
                                                     cluster.inside.de$gene <- rownames(cluster.inside.de)
                                                     return(cluster.inside.de)
                                                   }else{
                                                     return(list(NULL))
                                                   }
                                                 })

                      }, mc.cores=1)
                            names(de.list) <- unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME]))
                            de.list <- de.list[sapply(de.list, function(x){(!is.null(x[[1]]))})]
                            nm <- paste0(ident.use[1], ".vs.", ident.use[2])
                            return(de.list)
           }, mc.cores=1)

  names(all_de_list) <- sapply(lst, function(x) paste0(x[1], ".vs.", x[2]) )

  for (nm in names(all_de_list)) {
    logger.info("****processing go up&down %s", nm)
    de.list <- all_de_list[[nm]]
    go_ups <- get_go_up(de.list)
    all_goup_list[[nm]] <- go_ups
    go_downs <- get_go_down(de.list)
    all_godown_list[[nm]] <- go_downs
    dego_dump(nm, de.list, go_ups, go_downs)
    rm(de.list)
    rm(go_ups)
    rm(go_downs)
  }

  store_list <- list(all_de_list, all_goup_list, all_godown_list)
  names(store_list) <- c("de", "goup", "godown")

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("dego_stage_%s.Rds", DEFUALT_CLUSTER_NAME))
    save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("dego_stage_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("dego_stage_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }
  rm(store_list)
  rm(all_de_list)
  rm(all_goup_list)
  rm(all_godown_list)
  return(list(scrna, ret_code))
}
generate_scrna_dego_stage_vsRest <- function(scrna){
  ######  need to check if the DEFUALT_CLUSTER_NAME existed or not
  ret_code = 0
  lst <- unique(stage_lst)
  all_de_list = list()
  all_goup_list = list()
  all_godown_list = list()
  all_de_list <- mclapply(lst, function(ident.use){
                            logger.info("****processing %s vs Rest", ident.use)
                            Idents(scrna) <- DEFUALT_CLUSTER_NAME
                            de.list <- list()
                            de.list = mclapply(unique(sort(scrna@meta.data[, DEFUALT_CLUSTER_NAME])), function(clusters){
                                                 logger.info("processing cluster %s", as.character(clusters))
                                                 sub_idents = which(scrna@meta.data[, DEFUALT_CLUSTER_NAME] == clusters)
                                                 if(length(sub_idents) == 0){
                                                   return(list(NULL))
                                                   next
                                                 }
                                                 a_sub.subset <- subset(scrna, cells = sub_idents)
                                                 Idents(a_sub.subset) <- "stage"

                                                 cluster.inside.de <- list(NULL)
                                                 tryCatch({
                                                   cluster.inside.de <- FindMarkers(a_sub.subset, ident.1= ident.use, logfc.threshold=0)
                                                 }, error=function(cond) {
                                                   logger.warn(cond)
                                                   cluster.inside.de <- list(NULL)
                                                 },finally={
                                                   if (length(cluster.inside.de) != 0){
                                                     cluster.inside.de$gene <- rownames(cluster.inside.de)
                                                     return(cluster.inside.de)
                                                   }else{
                                                     return(list(NULL))
                                                   }
                                                 })

                      }, mc.cores=1)
                            names(de.list) <- unique(sort(scrna@meta.data[, DEFUALT_CLUSTER_NAME]))
                            de.list <- de.list[sapply(de.list, function(x){(!is.null(x[[1]]))})]
                            nm <- ident.use
                            return(de.list)
           }, mc.cores=1)

  save_object(
    object = all_de_list,
    file_name = file.path(SAVE_DIR,"all_de_list.Rds"),
    file_format = COMPRESSION_FORMAT
  )
  names(all_de_list) <- sapply(lst, function(x) paste0(x, ".vs.Rest") )

  for(nm in names(all_de_list)) {
    logger.info("****processing go up&down %s", nm)
    de.list <- all_de_list[[nm]]
    go_ups <- get_go_up(de.list)
    all_goup_list[[nm]] <- go_ups
    go_downs <- get_go_down(de.list)
    all_godown_list[[nm]] <- go_downs
    dego_dump(nm, de.list, go_ups, go_downs)
    rm(de.list)
    rm(go_ups)
    rm(go_downs)
  }

  store_list <- list(all_de_list, all_goup_list, all_godown_list)
  names(store_list) <- c("de", "goup", "godown")

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("dego_stage_vsRest_%s.Rds", DEFUALT_CLUSTER_NAME))
    save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("dego_stage_vsRest_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("dego_stage_vsRest_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }

  rm(store_list)
  rm(all_de_list)
  rm(all_goup_list)
  rm(all_godown_list)

  return(list(scrna, ret_code))
}


generate_scrna_go <- function(scrna){
  ret_code = 0

  partition = sprintf("de_%s", DEFUALT_CLUSTER_NAME)

  if(!ALLINONE){
    de.list <- seutools_partition(scrna, partition, SAVE_DIR, allinone=FALSE)
  }else{
    de.list <- scrna@tools[[partition]]
  }
  go.up.list <- get_go_up(de.list)
  go.down.list <- get_go_down(de.list)
  store_list <- list(go.up.list, go.down.list)
  names(store_list) <- c("goup", "godown")

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("go_%s.Rds", DEFUALT_CLUSTER_NAME))
    save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("go_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("go_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }

  golist_xls(go.up.list,  sprintf("goup_%s.xlsx", DEFUALT_CLUSTER_NAME))
  golist_xls(go.down.list, sprintf("godown_%s.xlsx", DEFUALT_CLUSTER_NAME))

  rm(go.up.list)
  rm(go.down.list)
  return(list(scrna, ret_code))
}

generate_scrna_MSigDB_geneset <- function(scrna){
    #http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/msigdb.v7.4.symbols.gmt
  ret_code = 0
  suppressPackageStartupMessages(require(GSEABase))
  tryCatch(
           {
             HS_Gmt <- getGmt(gzfile(MSigDB_GENESET_HUMAN_GMT_FILE))
             Gmt <- HS_Gmt
             if(SPECIES == "Mouse"){
               all_mappings <- read.csv(file="external/Human2Mouse_mappings.tsv", sep="\t")
               for(idx in 1:length(HS_Gmt)){
                  geneIds <- HS_Gmt[[idx]]@geneIds
                  mm_geneIds <-all_mappings[which(all_mappings$HGNC.symbol %in% geneIds), "MGI.symbol"]
                  Gmt[[idx]]@geneIds <- mm_geneIds
               }
             }
             assertthat::assert_that(all(MSigDB_Geneset_names %in% names(Gmt)))
             subGmt <- Gmt[MSigDB_Geneset_names]

             for(nm in names(subGmt)){
                 geneIds <- subGmt[[nm]]@geneIds
                 geneIds <- intersect(geneIds, rownames(scrna))
                 ## scrna@assays$RNA@data
                 scrna <- AddModuleScore(object = scrna, features = list(geneIds), name = nm, assay="RNA")
                 scrna@meta.data[, nm] <- scrna@meta.data[, paste0(nm, 1)]
                 scrna@meta.data[, paste0(nm, 1)] <- NULL
             }
             scrna@tools$genesets <- names(subGmt)
             rm(Gmt)
             rm(subGmt)
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

generate_scrna_progeny <- function(scrna){
  ret_code = 0
  assertthat::assert_that(SPECIES == "Human" | SPECIES == "Mouse")
  Idents(scrna) <- DEFUALT_CLUSTER_NAME
  scrna <- progeny::progeny(scrna, scale=FALSE, organism=SPECIES, top=500, perm=1,return_assay = TRUE)

  da <- DefaultAssay(scrna)
  DefaultAssay(scrna) <-'progeny'
  pws <- rownames(scrna@assays$progeny)
  res <- list()
  Idents(scrna) <- DEFUALT_CLUSTER_NAME
  if (is.factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- droplevels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
  }else{
    c_names <-unique(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME],
                                                      levels=sort(help_sort_func(c_names)))
  }

  for(i in levels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    g <- as.character(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    g[!(g==i)] <- "others"
    g <- factor(g, levels=c(i, "others"))
    res[[i]] = scran::findMarkers(as.matrix(scrna@assays$progeny@data), g)[[1]]
    res[[i]] <- as.data.frame(res[[i]])
    r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(scrna@assays$progeny@data[pw,]), g))
    nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
    names(r) <- nms
    res[[i]][nms, "r"] <- r
    res[[i]] <- res[[i]][nms, ]
  }

  for (cl in names(res)) {
    res[[cl]]$pathway <- rownames(res[[cl]])
    res[[cl]]$CellType <- cl
    colnames(res[[cl]]) <-  c("Top","p.value","FDR", "summary.logFC","logFC","r","pathway","CellType")
  }
  res_df <- do.call("rbind", res)
  res_df$tag <- sapply(res_df$FDR, function(pval) {
      if(pval< 0.001) {
      txt <- "***"
      } else if (pval < 0.01) {
      txt <- "**"
      } else if (pval < 0.05) {
      txt <- "*"
      } else {
      txt <- ""
      }
      return(txt)
  })
  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", glue("progeny_{DEFUALT_CLUSTER_NAME}.Rds"))
    save_object(res_df,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[glue("progeny_{DEFUALT_CLUSTER_NAME}")]] <- fname
  }else{
     scrna@tools[[glue("progeny_{DEFUALT_CLUSTER_NAME}")]] <- res_df
  }
  rm(res_df)

  DefaultAssay(scrna) <- da
  return(list(scrna, ret_code))
}

generate_scrna_progeny_stage <- function(scrna){

  ret_code = 0
  if("progeny" %ni% names(scrna@assays)){
    scrna <- progeny::progeny(scrna, scale=FALSE, organism=SPECIES, top=500, perm=1,return_assay = TRUE)
  }

  if (!ALLINONE){
    fname = file.path(SAVE_DIR, "assays", "progeny.Rds")
    # construct list
    assay_info <- list(
      name = "progeny",
      assay = scrna[["progeny"]],
      fname = fname,
      meta = scrna@meta.data,
      info = "progeny pathway scores")
    #save to disk
    save_object(assay_info, fname, file_format = COMPRESSION_FORMAT)
    if ("assay_info" %ni% names(scrna@tools)){
       scrna@tools[["assay_info"]] <- list()
    }
    scrna@tools[["assay_info"]][["progeny"]] <- fname
    ## remove from memory
    rm(assay_info)
    # small assay still keep it here
    #scrna[["progeny"]] <- NULL
    gc()
  }



  conds <- scrna@tools$meta_order$stage
  m <- combn(conds, 2)
  n = length(m)/2
  lst <- vector("list", n)
  for (i in 1:n){
    lst[[i]] <- m[1:2, i]
  }
  if (is.factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- droplevels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
  }else{
    c_names <-unique(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME],
                                                      levels=sort(help_sort_func(c_names)))
  }

  pws <- rownames(scrna@assays$progeny)
  vs_df_list <- list()
  for(apair in lst){
    vs1 <- apair[1]
    vs2 <- apair[2]

    Idents(scrna) <- 'stage'
    ###
    res <- list()
    for(i in levels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
        cells1 <- which(scrna@meta.data[, DEFUALT_CLUSTER_NAME]==i & (scrna@meta.data$stage %in% vs1))
        cells2 <- which(scrna@meta.data[, DEFUALT_CLUSTER_NAME]==i & (scrna@meta.data$stage %in% vs2))
        if( length(cells1) < 2 | length(cells2) < 2){
           next
        }
        a_sub = subset(scrna, cells=c(cells1, cells2))
        g <- as.character(a_sub@meta.data$stage)
        g <- factor(g, levels=c(vs1, vs2))
        res[[i]] = scran::findMarkers(as.matrix(a_sub@assays$progeny@data), g)[[1]]
        res[[i]] <- as.data.frame(res[[i]])
        r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$progeny@data[pw,]), g))
        nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
        names(r) <- nms
        res[[i]][nms, "r"] <- r
        res[[i]] <- res[[i]][nms, ]

    }

    for (cl in names(res)) {
        res[[cl]]$pathway <- rownames(res[[cl]])
        res[[cl]]$CellType <- cl
        colnames(res[[cl]]) <-  c("Top","p.value","FDR", "summary.logFC","logFC","r","pathway","CellType")

    }
    res_df <- do.call("rbind", res)
#    res_df$FDR <- p.adjust(res_df$p_val, method="fdr")
    res_df$tag <- sapply(res_df$FDR, function(pval) {
        if(pval< 0.001) {
        txt <- "***"
        } else if (pval < 0.01) {
        txt <- "**"
        } else if (pval < 0.05) {
        txt <- "*"
        } else {
        txt <- ""
        }
        return(txt)
    })

    vs_df_list[[glue("{vs1}.vs.{vs2}")]] <- res_df
    rm(res_df)

  }
  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", glue("progeny_stage_{DEFUALT_CLUSTER_NAME}.Rds"))
     save_object( vs_df_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[glue("progeny_stage_{DEFUALT_CLUSTER_NAME}")]] <- fname
  }else{
     scrna@tools[[glue("progeny_stage_{DEFUALT_CLUSTER_NAME}")]] <- vs_df_list
  }
  rm(vs_df_list)
  return(list(scrna, ret_code))
}


generate_scrna_pathway_stage <- function(scrna){
  ret_list <- get_pathway_comparison(scrna, slot="stage")
  scrna <- ret_list[[1]]
  ret_code <- ret_list[[2]]
  return(list(scrna, ret_code))
}

generate_scrna_pathway_name <- function(scrna){
  ret_list <- get_pathway_comparison(scrna, slot="name")
  scrna <- ret_list[[1]]
  ret_code <- ret_list[[2]]
  return(list(scrna, ret_code))
}

generate_scrna_pathway_stage_vsRest <- function(scrna){
  ret_list <- get_pathway_comparison(scrna, slot="stage_vsRest")
  scrna <- ret_list[[1]]
  ret_code <- ret_list[[2]]
  return(list(scrna, ret_code))
}

generate_scrna_kegg <- function(scrna){
  ret_code = 0
  if(!ALLINONE){
    de.list <- seutools_partition(scrna, sprintf("de_%s", DEFUALT_CLUSTER_NAME), SAVE_DIR, allinone=FALSE)
  }else{
    de.list <- scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]]
  }
  kegg.up.list <- get_kegg_up(de.list)
  kegg.down.list <- get_kegg_down(de.list)
  store_list <- list(kegg.up.list, kegg.down.list)
  names(store_list) <- c("keggup", "keggdown")
  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("kegg_%s.Rds", DEFUALT_CLUSTER_NAME))
     save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("kegg_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("kegg_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }

  golist_xls(kegg.up.list,  sprintf("keggup_%s.xlsx", DEFUALT_CLUSTER_NAME))
  golist_xls(kegg.down.list, sprintf("keggdown_%s.xlsx", DEFUALT_CLUSTER_NAME))
  rm(store_list)
  rm(kegg.up.list)
  rm(kegg.down.list)

  return(list(scrna, ret_code))
}


generate_scrna_reactome <- function(scrna){
  ret_code = 0
  if(!ALLINONE){
    de.list <- seutools_partition(scrna, sprintf("de_%s", DEFUALT_CLUSTER_NAME), SAVE_DIR, allinone=FALSE)
  }else{
    de.list <- scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]]
  }
  reactome.up.list <- get_reactome_up(de.list)
  reactome.down.list <- get_reactome_down(de.list)
  store_list <- list(reactome.up.list, reactome.down.list)
  names(store_list) <- c("reactomeup", "reactomedown")

  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("reactome_%s.Rds", DEFUALT_CLUSTER_NAME))
     save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("reactome_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("reactome_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }
  golist_xls(reactome.up.list,  sprintf("reactomeup_%s.xlsx", DEFUALT_CLUSTER_NAME))
  golist_xls(reactome.down.list, sprintf("reactomedown_%s.xlsx", DEFUALT_CLUSTER_NAME))

  rm(store_list)
  rm(reactome.up.list)
  rm(reactome.down.list)

  return(list(scrna, ret_code))
}


generate_scrna_hallmark <- function(scrna){
  ret_code = 0
  if(!ALLINONE){
    de.list <- seutools_partition(scrna, sprintf("de_%s", DEFUALT_CLUSTER_NAME), SAVE_DIR, allinone=FALSE)
  }else{
    de.list <- scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]]
  }
  hallmark.up.list <- get_hallmark_up(de.list)
  hallmark.down.list <- get_hallmark_down(de.list)
  store_list <- list(hallmark.up.list, hallmark.down.list)
  names(store_list) <- c("hallmarkup", "hallmarkdown")


  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", sprintf("hallmark_%s.Rds", DEFUALT_CLUSTER_NAME))
     save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[sprintf("hallmark_%s", DEFUALT_CLUSTER_NAME)]] <- fname
  }else{
     scrna@tools[[sprintf("hallmark_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
  }

  golist_xls(hallmark.up.list,  sprintf("hallmarkup_%s.xlsx", DEFUALT_CLUSTER_NAME))
  golist_xls(hallmark.down.list, sprintf("hallmarkdown_%s.xlsx", DEFUALT_CLUSTER_NAME))
  rm(store_list)
  rm(hallmark.up.list)
  rm(hallmark.down.list)

  return(list(scrna, ret_code))
}


get_go_up <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  go.up.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  if(orgdb == "NotSupport"){
    logger.error("not support species for GO!!!")
    stop("not support species for GO!!!")
  }
  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id,  date()))

    genes_e = NULL
    genes.sorted = NULL
    tryCatch({
      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      go.up.list[(id)] = list(NULL)
      print("null")
      next
    }

    go <- enrichGO(genes_e$ENTREZID, OrgDb = orgdb, ont='ALL',pAdjustMethod = 'BH',
                   pvalueCutoff = 1, qvalueCutoff = 1,keyType = 'ENTREZID', readable=TRUE)
    #barplot(go,showCategory=20, drop=T, title = paste("GO analysis for LSK UP genes for cluster", id, sep = " "))
    go.up.list[[id]] = go
  }
  if(length(go.up.list) == 0){
    return(list())
  }
  go.up.list <- go.up.list[sapply(go.up.list, function(x){!is.null(x)})]
  return(go.up.list)
}


get_go_down <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  go.down.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  if(orgdb == "NotSupport"){
    logger.error("not support species for GO!!!")
    stop("not support species for GO!!!")
  }
  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id, date()))
    genes_e = NULL
    genes.sorted = NULL

    tryCatch({

      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      go.down.list[(id)] = list(NULL)
      print("null")
      next
    }


    go <- enrichGO(genes_e$ENTREZID, OrgDb = orgdb, ont='ALL',pAdjustMethod = 'BH',
                   pvalueCutoff = 1, qvalueCutoff = 1,keyType = 'ENTREZID', readable=TRUE)
    #barplot(go,showCategory=20, drop=T, title = paste("GO analysis for LSK Down genes for cluster", id, sep = " "))
    go.down.list[[id]] = go
  }
  if(length(go.down.list) == 0){
    return(list())
  }
  go.down.list <- go.down.list[sapply(go.down.list, function(x){!is.null(x)})]
  return(go.down.list)
}

get_kegg_up <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  kegg.up.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  keggorgan <- switch(SPECIES, "Mouse"="mmu",
                      "Human"="hsa",
                      "NotSupport")

  if(orgdb == "NotSupport" | keggorgan == "NotSupport"){
    logger.error("not support species for KEGG!!!")
    stop("not support species for KEGG!!!")
  }

  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id,  date()))

    genes_e = NULL
    genes.sorted = NULL
    tryCatch({
      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      kegg.up.list[(id)] = list(NULL)
      print("null")
      next
    }

    kegg <- enrichKEGG(gene= genes_e$ENTREZID,
                       organism     = keggorgan,
                       pvalueCutoff = 1)

    if(is.null(kegg)){
      kegg.up.list[(id)] = list(NULL)
      print("null")
      next
    }
    kegg <- setReadable(kegg, orgdb, keyType = "ENTREZID")
    kegg.up.list[[id]] = kegg
  }
  if(length(kegg.up.list) == 0){
    return(list())
  }
  kegg.up.list <- kegg.up.list[sapply(kegg.up.list, function(x){!is.null(x)})]
  return(kegg.up.list)
}


get_kegg_down <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  kegg.down.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  keggorgan <- switch(SPECIES, "Mouse"="mmu",
                      "Human"="hsa",
                      "NotSupport")

  if(orgdb == "NotSupport" | keggorgan == "NotSupport"){
    logger.error("not support species for KEGG!!!")
    stop("not support species for KEGG!!!")
  }
  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id, date()))
    genes_e = NULL
    genes.sorted = NULL

    tryCatch({

      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      kegg.down.list[(id)] = list(NULL)
      print("null")
      next
    }

    kegg <- enrichKEGG(gene= genes_e$ENTREZID,
                       organism     = keggorgan,
                       pvalueCutoff = 1)

    if(is.null(kegg)){
      kegg.down.list[(id)] = list(NULL)
      print("null")
      next
    }

    kegg <- setReadable(kegg, orgdb, keyType = "ENTREZID")
    kegg.down.list[[id]] = kegg
  }
  if(length(kegg.down.list) == 0){
    return(list())
  }
  kegg.down.list <- kegg.down.list[sapply(kegg.down.list, function(x){!is.null(x)})]
  return(kegg.down.list)
}

get_reactome_up <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(ReactomePA))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  reactome.up.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  reactomeorgan <- switch(SPECIES, "Mouse"="mouse",
                          "Human"="human",
                          "NotSupport")

  if(orgdb == "NotSupport" | reactomeorgan == "NotSupport"){
    logger.error("not support species for reactome!!!")
    stop("not support species for reactome!!!")
  }

  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id,  date()))

    genes_e = NULL
    genes.sorted = NULL
    tryCatch({
      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      reactome.up.list[(id)] = list(NULL)
      print("null")
      next
    }
    reactome <- enrichPathway(gene=genes_e$ENTREZID, organism = reactomeorgan, readable=T)
    if(is.null(reactome)){
      reactome.up.list[(id)] = list(NULL)
    }else{
      reactome.up.list[[id]] = reactome
    }
  }
  if(length(reactome.up.list) == 0){
    return(list())
  }
  reactome.up.list <- reactome.up.list[sapply(reactome.up.list, function(x){!is.null(x)})]
  return(reactome.up.list)
}


get_reactome_down <- function(de.list){
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(ReactomePA))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  reactome.down.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  reactomeorgan <- switch(SPECIES, "Mouse"="mouse",
                          "Human"="human",
                          "NotSupport")

  if(orgdb == "NotSupport" | reactomeorgan == "NotSupport"){
    logger.error("not support species for reactome!!!")
    stop("not support species for reactome!!!")
  }
  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id, date()))
    genes_e = NULL
    genes.sorted = NULL

    tryCatch({

      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      reactome.down.list[(id)] = list(NULL)
      print("null")
      next
    }

    reactome <- enrichPathway(gene=genes_e$ENTREZID, organism = reactomeorgan, readable=T)
    if(is.null(reactome)){
      reactome.down.list[(id)] = list(NULL)
    }else{
      reactome.down.list[[id]] = reactome
    }
    reactome.down.list[[id]] = reactome
  }
  if(length(reactome.down.list) == 0){
    return(list())
  }
  reactome.down.list <- reactome.down.list[sapply(reactome.down.list, function(x){!is.null(x)})]
  return(reactome.down.list)
}


get_hallmark_up <- function(de.list){
  suppressPackageStartupMessages(require(msigdbr))
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  hallmark.up.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")

  hallmarkorgan <- switch(SPECIES, "Mouse"="Mus musculus",
                          "Human"="Homo sapiens",
                          "NotSupport")

  if(orgdb == "NotSupport" | hallmarkorgan == "NotSupport"){
    logger.error("not support species for hallmark!!!")
    stop("not support species for hallmark!!!")
  }

  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC > 0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id,  date()))

    genes_e = NULL
    genes.sorted = NULL
    tryCatch({
      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      hallmark.up.list[(id)] = list(NULL)
      print("null")
      next
    }

    m_t2g <- msigdbr(species = hallmarkorgan, category = "H") %>%
      dplyr::select(gs_name, entrez_gene)

    if(length(intersect(genes_e$ENTREZID,  m_t2g$entrez_gene)) == 0){
      hallmark.up.list[(id)] = list(NULL)
      print("null")
      next
    }

    hallmark <- enricher(genes_e$ENTREZID, TERM2GENE=m_t2g, pvalueCutoff=1)
    hallmark <- setReadable(hallmark, orgdb, keyType = "ENTREZID")
    hallmark.up.list[[id]] = hallmark
  }
  if(length(hallmark.up.list) == 0){
    return(list())
  }
  hallmark.up.list <- hallmark.up.list[sapply(hallmark.up.list, function(x){!is.null(x)})]
  return(hallmark.up.list)
}


get_hallmark_down <- function(de.list){
  suppressPackageStartupMessages(require(msigdbr))
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  hallmark.down.list = list()
  help_sort_func <- ifelse(all.is.numeric(names(de.list)), as.numeric, function(x){x})
  orgdb <- switch(SPECIES, "Mouse"="org.Mm.eg.db",
                  "Human"="org.Hs.eg.db",
                  "NotSupport")
  hallmarkorgan <- switch(SPECIES, "Mouse"="Mus musculus",
                          "Human"="Homo sapiens",
                          "NotSupport")

  if(orgdb == "NotSupport" | hallmarkorgan == "NotSupport"){
    logger.error("not support species for hallmark!!!")
    stop("not support species for hallmark!!!")
  }
  for (id in sort(help_sort_func(names(de.list)))) {
    id <- as.character(id)
    if("avg_logFC" %in% names(de.list[[id]])){ ## compatible for seurat3
        de.list[[id]]$avg_log2FC <- de.list[[id]]$avg_logFC/log(2)
    }
    genes = de.list[[id]]$gene[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_log2FC < -0.36 & de.list[[id]]$p_val_adj < 0.05]
    logger.info(paste("cluster: ", id, date()))
    genes_e = NULL
    genes.sorted = NULL

    tryCatch({
      genes.sorted <- genes[order(pvals)][1:min(100, length(genes))]
      genes_e <- bitr(genes.sorted, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=orgdb);
    }, error = function(cond) {
      logger.warn(cond)
    })

    if(is.null(genes_e)){
      hallmark.down.list[(id)] = list(NULL)
      print("null")
      next
    }

    m_t2g <- msigdbr(species = hallmarkorgan, category = "H") %>%
      dplyr::select(gs_name, entrez_gene)

    if(length(intersect(genes_e$ENTREZID,  m_t2g$entrez_gene)) == 0){
      hallmark.down.list[(id)] = list(NULL)
      print("null")
      next
    }

    hallmark <- enricher(genes_e$ENTREZID, TERM2GENE=m_t2g, pvalueCutoff=1)
    hallmark <- setReadable(hallmark, orgdb, keyType = "ENTREZID")
    hallmark.down.list[[id]] = hallmark
  }
  if(length(hallmark.down.list) == 0){
    return(list())
  }
  hallmark.down.list <- hallmark.down.list[sapply(hallmark.down.list, function(x){!is.null(x)})]
  return(hallmark.down.list)
}


get_pathway_comparison <- function(scrna, slot){
  ret_code <- 0

  dego_name <- glue("dego_{slot}_{DEFUALT_CLUSTER_NAME}")
  if (!(dego_name %in% names(scrna@tools) )){
    logger.error(glue("{dego_name} not hasn't been calculated!!!"))
    ret_code <- -1
    return(list(scrna, ret_code))
  }

  all_keggup_list = list()
  all_keggdown_list = list()
  all_hallmarkup_list = list()
  all_hallmarkdown_list = list()
  all_reactomeup_list = list()
  all_reactomedown_list = list()

  if(!ALLINONE){
    all_dego_list  <- seutools_partition(scrna, dego_name, SAVE_DIR, allinone=FALSE)
  }else{
    all_dego_list <- scrna@tools[[dego_name]]
  }
  all_de_list <- all_dego_list[["de"]]
  for(nm in names(all_de_list)) {
    logger.info("****processing pathway up&down %s", nm)
    de.list <- all_de_list[[nm]]

    kegg_ups <- get_kegg_up(de.list)
    all_keggup_list[[nm]] <- kegg_ups
    kegg_downs <- get_kegg_down(de.list)
    all_keggdown_list[[nm]] <- kegg_downs

    hallmark_ups <- get_hallmark_up(de.list)
    all_hallmarkup_list[[nm]] <- hallmark_ups
    hallmark_downs <- get_hallmark_down(de.list)
    all_hallmarkdown_list[[nm]] <- hallmark_downs


    reactome_ups <- get_reactome_up(de.list)
    all_reactomeup_list[[nm]] <- reactome_ups
    reactome_downs <- get_reactome_down(de.list)
    all_reactomedown_list[[nm]] <- reactome_downs

    pathway_dump(nm, "KEGG", kegg_ups, kegg_downs)
    pathway_dump(nm, "hallmark", hallmark_ups, hallmark_downs)
    pathway_dump(nm, "Reactome", reactome_ups, reactome_downs)

    rm(kegg_ups)
    rm(kegg_downs)
    rm(hallmark_ups)
    rm(hallmark_downs)
    rm(reactome_ups)
    rm(reactome_downs)
  }
  store_list <- list(all_keggup_list, all_keggdown_list,
                     all_hallmarkup_list, all_hallmarkdown_list,
                     all_reactomeup_list, all_reactomedown_list)
  names(store_list) <- c("keggup", "keggdown",
                         "hallmarkup", "hallmarkdown",
                         "reactomeup", "reactomedown")


  if(!ALLINONE){
    fname = file.path(SAVE_DIR, "partition", glue("pathway_{slot}_{DEFUALT_CLUSTER_NAME}.Rds"))
     save_object(store_list,
                 file_name = fname,
                 file_format = COMPRESSION_FORMAT)
     scrna@tools[[glue("pathway_{slot}_{DEFUALT_CLUSTER_NAME}")]] <- fname
  }else{
     scrna@tools[[glue("pathway_{slot}_{DEFUALT_CLUSTER_NAME}")]] <- store_list
  }
  rm(store_list)
  rm(all_keggup_list, all_keggdown_list,
     all_hallmarkup_list, all_hallmarkdown_list,
     all_reactomeup_list, all_reactomedown_list)


  return(list(scrna, ret_code))
}

# Sanity Function
# This function is called every time a Seurat object is loaded or conf_main processes a new step.

# We have the singleton slots. When we run the pipeline for a single sample and use the comparing features, we need to check for these slots.
# But we do not want to check for these slots if we use the comparing features for a set of samples, i.e. in the normal use.
# So, if conf does not contain scrna_phase_singleton, scrna_sltn_batch_clustering or scrna_singleton_clustering, we remove the singleton slots from the attributes
sanity_singles <- c("scrna_phase_singleton|scrna_sltn_batch_clustering|scrna_singleton_clustering")
sanity_singles <- any(grepl(sanity_singles, names(conf)))
if(sanity_singles == FALSE){
  sanity_attributes <- sanity_attributes[!grepl("SINGLE", sanity_attributes[,1]),]
} else{
  # We have phase_singleton or one of the singleton steps.
  # So,	we remove the INTE_PCA and INTE_UMAP.
  sanity_attributes <- sanity_attributes[!grepl("INTE|DEFAULT_UMAP", sanity_attributes[,1]),]
}

sanity_function <- function(scrna, key){
  print(paste0("Checking slots for ", key))
  sanity_type <- key
  #sanity_slots_to_check <- grepl(sanity_type, sanity_attributes$type) | grepl("all", sanity_attributes$type)
  #sanity_attributes_use <- sanity_attributes #sanity_attributes[sanity_slots_to_check,]

  # We apply the check_attributes function below
  sanity_attributes_result <- lapply(sanity_attributes[,1], check_attributes, scrna = scrna)
  # We bind the results
  sanity_attributes_result <- do.call("rbind", sanity_attributes_result)
  scrna@misc[["sanity"]] <- sanity_attributes_result
  colnames(scrna@misc$sanity) <- key
  #scrna@misc[[paste0("sanity_", key)]] <- sanity_attributes_result
  print("Checked slot.")
  return(scrna)
}
check_attributes <- function(scrna, sanity_attribute){
  result <- tryCatch(
  {
    sanity_test <- eval(parse(text = sanity_attribute))
    sanity_result <- paste0("Present: ", sanity_attribute)
  },
  error = function(cond){
    sanity_result <- paste0("Missing: ", sanity_attribute)
  },
  finally = {
    #pass
  })
}

# DoubletFinder
generate_scrna_doublet_proportions <- function(scrna){
  suppressPackageStartupMessages(library(DoubletFinder))
  ret_code <- 1
  if(doublet_switch %in% c("on", "display")){
    doublet_formation_rate <- determine_doublet_proportions(scrna = scrna, doublet_lst = doublet_lst)
    #names(doublet_formation_rate) <- names(doublet_lst)
    scrna_lst <- SplitObject(scrna, split.by = "name")
    scrna_lst <- lapply(scrna_lst, function(x){
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      x <- ScaleData(x)
      x <- RunPCA(x)
      x <- RunUMAP(x, dims = 1:20)
    })
    bcmvn_lst <- lapply(scrna_lst, mc_pK_identification)
    pK_optimal <- lapply(X = bcmvn_lst, FUN = function(x){as.numeric(as.character(x[x$BCmetric == max(x$BCmetric),2]))})
    est_expected <- lapply(X = names(bcmvn_lst), FUN = mc_est_expected, scrnas = scrna_lst, doublet_rate = doublet_formation_rate)
    est_expected <- unlist(est_expected)

    scrna_list_doublets <- lapply(
      X = names(bcmvn_lst),
      FUN = mc_doubletFinder_v3,
      seurats = scrna_lst,
      pN = 0.25,
      pK_optimal = pK_optimal,
      est_expected = est_expected
    )

    names(scrna_list_doublets) <- names(bcmvn_lst)
    classifications <- c()
    pANN <- c()
    cells <- c()
    for(i in 1:length(scrna_list_doublets)){
      cells               <- c(names(scrna_list_doublets[[i]]$classifications), cells)
      classifications     <- c(as.character(scrna_list_doublets[[i]]$classifications), classifications)
      pANN                <- c(pANN, scrna_list_doublets[[i]]@meta.data[,paste0("pANN_0.25_", pK_optimal[i], "_", est_expected[i])])
    }
    classifications <- factor(classifications, levels = c("Singlet", "Doublet"))
    names(classifications) <- cells
    names(pANN) <- cells
    scrna <- AddMetaData(scrna, classifications, col.name = "Doublet_classifications")
    scrna <- AddMetaData(scrna, pANN, col.name = "pANN")

    if(length(unique(scrna$name)) > 1){
      data.list <- SplitObject(scrna, split.by = "name")
      data.list <- lapply(X = data.list, FUN = function(x) {
                          x <- NormalizeData(x)
                          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
      k.filter <- min(table(scrna$name))
      k.filter <- ifelse(k.filter < 200, k.filter, 200)
      anchors <- FindIntegrationAnchors(object.list = data.list, dims = INTEGRATED_DIM, scale=TRUE,
                                        k.filter = k.filter)## THIS IS CCA DIMENSIONS
      scrna_save <- IntegrateData(anchorset = anchors, dims = INTEGRATED_DIM, k.weight = k.filter) ## THIS IS PCA DIMENSION
      scrna_save <- ScaleData(scrna_save, verbose = FALSE)
      scrna_save <- RunPCA(scrna_save, npcs = 30, verbose = FALSE, reduction.name="DOUBLET_PCA")
      scrna_save <- RunUMAP(scrna_save, reduction = "DOUBLET_PCA", dims = 1:20, reduction.name="DOUBLET_UMAP")
    } else if(length(unique(scrna$name)) == 1){
      scrna_save <- ScaleData(scrna, verbose = FALSE)
      scrna_save <- FindVariableFeatures(scrna_save)
      scrna_save <- RunPCA(scrna_save, npcs = 30, verbose = FALSE, reduction.name="DOUBLET_PCA")
      scrna_save <- RunUMAP(scrna_save, reduction = "DOUBLET_PCA", dims = 1:20, reduction.name="DOUBLET_UMAP")
    }
    save_object(
      object = scrna_save,
      file_name = file.path(SAVE_DIR, "scrna_DoubletAnnotated.Rds"),
      file_format = COMPRESSION_FORMAT
    )
  }
  if(doublet_switch == "on"){
    scrna <- subset(scrna, Doublet_classifications == "Singlet")
  }
  ret_code <- 0
  print("Finished DoubletFinder")
  return(list(scrna, ret_code))
}

if(!interactive()) {
  conf_main()
}
