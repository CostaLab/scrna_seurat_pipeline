#!/usr/bin/Rscript

###Set VERSION
VERSION = "1.0.1"


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
    message("scRNA seurat pipeline\nVersion: \n\t", VERSION, "\n")
    quit(status = 0)
}


suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger




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

    parser <- add_option(parser, c("-d", "--dims4Integrate"), type="character", default="1:20",
                    help="Dims to keep for integrating [default %default]",
                    metavar="VECTOR")

    parser <- add_option(parser, c("-x", "--Dims4FindNeighbors"), type="character", default="1:12",
                    help="Dims kept for findNeighbors [default %default]",
                    metavar="VECTOR")

    parser <- add_option(parser, c("-r", "--clusterresolution"), type="numeric", default=0.5,
                    help="Resolution for clustering [default %default]",
                    metavar="numeric")

    parser <- add_option(parser, c("-a", "--defaultclustername"), type="character", default="seurat_clusters",
                    help="downstream analysis default cluster name  [default %default]",
                    metavar="character")


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
FINDNEIGHBORS_DIM     = eval(parse(text=pa$Dims4FindNeighbors))
CLUSTER_RESOLUTION    = pa$clusterresolution 
DEFUALT_CLUSTER_NAME  = pa$defaultclustername
CM_FORMAT             = pa$countmatrixformat
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

##--------------keep variables----------------------
pa$project_name       = PROJECT
pa$organ              = ORGAN 
pa$species            = SPECIES 
pa$mca_name           = MCA_NAME
pa$external_file      = ANNOTATION_EXTERNAL_FILE
pa$data_src           = data_src
pa$stage_lst          = stage_lst 


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



## executing plan
conf = conf[conf > 0]

if(conf[1] != 2 & names(conf)[1] != "scrna_rawdata"){
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
suppressPackageStartupMessages(library(dplyr))
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
registerDoParallel(cores=WORKER_NUM)


# Set up future for parallelization 
plan("multiprocess", workers = WORKER_NUM)
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
               scrna <- readRDS(file=file.path(SAVE_DIR, NRds))
               scrna@tools$parameter[[cur_time]] <- unlist(pa)
               scrna@tools$execution[[cur_time]] <- conf 
               logger.info(paste("finished loading", NRds))
           }else if(val==1){
               f_name = paste("generate_", key, sep="")
               f_call = paste(f_name, "(scrna)", sep="")
               logger.info(paste("executing", f_name))
               if(startsWith(key, "scrna")){
                   ret_list <-  eval(parse(text=f_call))
                   scrna  <- ret_list[[1]]
                   ret_code <- ret_list[[2]]
                   stopifnot(ret_code == 0)
                   if(key == "scrna_rawdata"){
                        scrna@tools$parameter[[cur_time]] <- unlist(pa)
                        scrna@tools$execution[[cur_time]] <- conf
                   }  
                   saveRDS(scrna, file=file.path(SAVE_DIR, NRds)) 
                   logger.info(paste("finished", f_name))
               }
           }else{
               logger.error("Wrong setting in config file in section RUN PARAMETERS\n Only 0 1 2 permitted")
               logger.error(traceback())
               stop("Exit 1")
           }
       }
       logger.info("===============Finished===============")
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
         SheetNames = names(flist))
       
     golist_xls(go_ups, sprintf("%s.goup_%s.xlsx",file_predix, DEFUALT_CLUSTER_NAME))
     golist_xls(go_downs,sprintf("%s.godown_%s.xlsx",file_predix, DEFUALT_CLUSTER_NAME))
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
    go.xls.list = lapply(lst, fortify, showCategory=Inf)
    names(go.xls.list) = names(go.list)
    WriteXLS(
        go.xls.list,
        file.path(CHARTS_DIR, fxls),
        SheetNames = names(go.xls.list))
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

	    scrna <- merge(x = data.list[[1]], y = unlist(data.list[2:length(data.list)]),
                       merge.data = FALSE, add.cell.ids = names(data_src), project = PROJECT)
        scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-|^MT-")
        scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rpl|^Rps|^RPL|^RPS")

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

         scrna <- ScaleData(scrna)
         scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes), reduction.name="BCELLCYCLE_PCA")
         scrna <- CellCycleScoring(scrna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

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
		scrna <- ScaleData(scrna, vars.to.regress = c("CC.Difference"), features = rownames(scrna))
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
generate_scrna_excludeRegressOut<- function(scrna){
    ret_code = 0
	tryCatch(
	{
		scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA","percent.exclude"), features = rownames(scrna))
    	scrna <- RunPCA(scrna, features = VariableFeatures(scrna), nfeatures.print = 10, reduction.name="EXCLUDE_PCA")
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

generate_scrna_RegressOutAll <- function(scrna){
    ret_code = 0
	tryCatch(
	{
        ### Scale data add exclude
        scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
		scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA", "percent.exclude","CC.Difference"), features = rownames(scrna))
    	scrna <- RunPCA(scrna, features = VariableFeatures(scrna), nfeatures.print = 10, reduction.name="RegressAll_PCA")
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

generate_scrna_integration <- function(scrna){
    ret_code = 0
	tryCatch(
	{
	      bak_tools <- scrna@tools
        data.list <- SplitObject(scrna, split.by = "name")
        data.list <- lapply(X = data.list, FUN = function(x) {
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

		    anchors <- FindIntegrationAnchors(object.list = data.list, dims = INTEGRATED_DIM)## THIS IS CCA DIMENSIONS 
    	  scrna <- IntegrateData(anchorset = anchors, dims = INTEGRATED_DIM) ## THIS IS PCA DIMENSION
        ## keep the order of name and stage
        scrna$name <- factor(scrna$name, levels=names(data_src))
        scrna$stage <- factor(scrna$stage, levels=unique(stage_lst[names(data_src)]))
        scrna@tools$parameter <- bak_tools$parameter
        scrna@tools$execution <- bak_tools$execution

	},
	error=function(cond) {
        ret_code <<- -1
		logger.error(cond)
        logger.error(traceback())
	},

	finally={    
        rm(anchors)
        rm(data.list)
		return(list(scrna, ret_code))
	}) 
    return(list(scrna, ret_code))
}



generate_scrna_ScaleIntegration <- function(scrna){
    ret_code = 0
	tryCatch(
	{
		DefaultAssay(scrna) <- "integrated"
    	# Run the standard workflow for visualization and clustering
    	scrna <- ScaleData(scrna, verbose = FALSE)
    	scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE, reduction.name="INTE_PCA")

        scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
    	scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE, reduction.name="INTE_PCA_UMI")


    	scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA", "percent.exclude"), verbose = FALSE)
    	scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE, reduction.name="INTE_PCA_EXCLUDE")

        scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
    	scrna <- ScaleData(scrna, vars.to.regress = c("nCount_RNA", "percent.exclude","CC.Difference"), verbose = FALSE)
    	scrna <- RunPCA(scrna, npcs = 30, verbose = FALSE, reduction.name="INTE_PCA_ALL")

    	scrna <- RunUMAP(scrna, reduction = "INTE_PCA", dims = 1:20, reduction.name="INTE_UMAP")
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
        a_meta <- sprintf("integrated_snn_res.%s", CLUSTER_RESOLUTION)
        if(a_meta %in% colnames(scrna@meta.data)){
            scrna$seurat_clusters <- scrna@meta.data[, a_meta]
        }else{
		    scrna <- FindNeighbors(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM)
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

generate_scrna_batchclustering <- function(scrna){
    ret_code = 0
	tryCatch(
	{   
		scrna <- FindNeighbors(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM)
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

generate_scrna_fishertest_clusters <- function(scrna){
    ret_code = 0
    CLUSTER_TO_TEST <- DEFUALT_CLUSTER_NAME 
    STAGE_TO_TEST <- "stage"
    stages <- names(table(scrna@meta.data[, STAGE_TO_TEST]))
    
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
    for (x in stages){
         
         new_column_nm <- paste0(x, "_excluded")
         df[, new_column_nm] <- rowSums(count.matrix) - df[, x]
    }

    for (x in stages){
        exclude_nm <- paste0(x, "_excluded")
        df[, paste0("o.", x)] <- sum(df[, x]) - df[, x]
        df[, paste0("o.", exclude_nm)] <- sum(df[, exclude_nm]) - df[, exclude_nm]

    }

    for (x in stages){
        exclude_nm <- paste0(x, "_excluded")
        o_x <- paste0("o.", x)
        o_e <- paste0("o.", exclude_nm)
        vec <- c(x, exclude_nm, o_x, o_e)
        nm_mtx <-  matrix(vec, nrow=2)
        for(cluster in df$Cluster){
            mtx <- matrix(unlist(df[cluster, vec]), nrow=2)
            ft <- fisher.test(mtx,workspace=1e9)
            df[cluster, paste0("pval_", x)] <- ft$p.value
            df[cluster, paste0("odds.ratio_", x)] <- ft$estimate
        }
        df[, paste0("pval.adjust_", x)] <- p.adjust(df[, paste0("pval_", x)],  method = "bonferroni")
    }
    scrna@tools[[sprintf("fishertest_%s", DEFUALT_CLUSTER_NAME)]] <- df
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
     return(list(scrna, ret_code)) 
}


generate_scrna_remove_recluster <- function(scrna){
     ret_code = 0
     Idents(scrna) <- "seurat_clusters"
     keeps <- setdiff(unique(scrna$seurat_clusters), scrna_remove_clusters)
     scrna <- subset(scrna, idents= keeps)
     ret_list <- generate_scrna_integration(scrna)
     scrna <- ret_list[[1]]
     ret_code <- ret_list[[2]]
     if(ret_code != 0){
        return(list(scrna, ret_code)) 
     }

     ret_list <- generate_scrna_ScaleIntegration(scrna)
     scrna <- ret_list[[1]]
     ret_code <- ret_list[[2]]
     if(ret_code != 0){
        return(list(scrna, ret_code)) 
     }

     scrna <- FindNeighbors(scrna, reduction = "INTE_PCA", dims = FINDNEIGHBORS_DIM)
     scrna <- FindClusters(scrna, resolution = CLUSTER_RESOLUTION_RANGE) ## 
     scrna$remove_recluster <- scrna$seurat_clusters
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
       sdf <- df[df$Mouse.Gene %in% names(score), c(glue("{SPECIES}.Gene"), "Cell.Type")] 
       sdf$score <- score[sdf[,glue("{SPECIES}.Gene")]] 
       itype <- aggregate(sdf$score, by=list(CellType=sdf$Cell.Type), 
                          FUN=function(x){return(mean(x)*log(length(x)+1))})
       celltype <- itype[which.max(itype$x), ]$CellType  
       return(celltype)
   }
   celltypes <- unlist(celltype.list)
   names(celltypes)<-colnames(mtx) 

   scrna$external_annotation <- unlist(celltypes)
   return(list(scrna, ret_code))

}


generate_scrna_markergenes <- function(scrna){
    ret_code = 0
	tryCatch(
	{
		DefaultAssay(scrna) <- "RNA"
        Idents(scrna) <- DEFUALT_CLUSTER_NAME 
    	de.df = FindAllMarkers(scrna) 

        cluster.de <- split(de.df, de.df$cluster)
        flist <- lapply(cluster.de, subset, subset = p_val_adj < 0.05)
        flist <-  flist[sapply(flist, function(m) nrow(m) >0)]
        WriteXLS(
                flist,
                file.path(CHARTS_DIR, sprintf("de_%s.xlsx", DEFUALT_CLUSTER_NAME)),
                SheetNames = names(flist))
        scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]] <- flist
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
        cluster.de.list <- foreach(i=CLUSTER_RESOLUTION_RANGE) %dopar%{
                    DefaultAssay(scrna) <- "RNA"
                    cluster_name <- sprintf("integrated_snn_res.%s", i)    	    
                    Idents(scrna) <- cluster_name 
                    de.df = FindAllMarkers(scrna)
                    cluster.de <- split(de.df, de.df$cluster)
        }
         
        scrna@tools[["de_batch"]] <- cluster.de.list
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
              cluster.inside.de <- FindMarkers(a_sub.subset, ident.1= ident.use[1])
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

            
      }, mc.cores=WORKER_NUM)
      names(de.list) <- unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME]))
      de.list <- de.list[sapply(de.list, function(x){(!is.null(x[[1]]))})]
      nm <- paste0(ident.use[1], ".vs.", ident.use[2])
      return(de.list)
    }, mc.cores=WORKER_NUM)
    
    names(all_de_list) <- sapply(lst, function(x) paste0(x[1], ".vs.", x[2]) )
  
    foreach(nm = names(all_de_list)) %do% {
      logger.info("****processing go up&down %s", nm)
      de.list <- all_de_list[[nm]]
      go_ups <- get_go_up(de.list)
      all_goup_list[[nm]] <- go_ups
      go_downs <- get_go_down(de.list)
      all_godown_list[[nm]] <- go_downs
      dego_dump(nm, de.list, go_ups, go_downs)
    }

    store_list <- list(all_de_list, all_goup_list, all_godown_list)
    names(store_list) <- c("de", "goup", "godown")
    scrna@tools[[sprintf("dego_name_%s", DEFUALT_CLUSTER_NAME)]] <- store_list 

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
              cluster.inside.de <- FindMarkers(a_sub.subset, ident.1= ident.use[1])
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
         
      }, mc.cores=WORKER_NUM)
      names(de.list) <- unique(sort(a_sub@meta.data[, DEFUALT_CLUSTER_NAME]))
      de.list <- de.list[sapply(de.list, function(x){(!is.null(x[[1]]))})]
      nm <- paste0(ident.use[1], ".vs.", ident.use[2])
      return(de.list)
    }, mc.cores=WORKER_NUM)

    names(all_de_list) <- sapply(lst, function(x) paste0(x[1], ".vs.", x[2]) )

    foreach(nm = names(all_de_list)) %do% {
      logger.info("****processing go up&down %s", nm)
      de.list <- all_de_list[[nm]]
      go_ups <- get_go_up(de.list)
      all_goup_list[[nm]] <- go_ups
      go_downs <- get_go_down(de.list)
      all_godown_list[[nm]] <- go_downs
      dego_dump(nm, de.list, go_ups, go_downs)
    }

   store_list <- list(all_de_list, all_goup_list, all_godown_list)
   names(store_list) <- c("de", "goup", "godown")
   scrna@tools[[sprintf("dego_stage_%s", DEFUALT_CLUSTER_NAME)]] <- store_list 

   return(list(scrna, ret_code))
}

generate_scrna_go <- function(scrna){
    ret_code = 0
    de.list <- scrna@tools[[sprintf("de_%s", DEFUALT_CLUSTER_NAME)]]
    go.up.list <- get_go_up(de.list)
    go.down.list <- get_go_down(de.list)
    store_list <- list(go.up.list, go.down.list)
    names(store_list) <- c("goup", "godown")
    scrna@tools[[sprintf("go_%s", DEFUALT_CLUSTER_NAME)]] <- store_list
    golist_xls(go.up.list,  sprintf("goup_%s.xlsx", DEFUALT_CLUSTER_NAME))
    golist_xls(go.down.list, sprintf("godown_%s.xlsx", DEFUALT_CLUSTER_NAME))

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
        genes = de.list[[id]]$gene[de.list[[id]]$avg_logFC > 0 & de.list[[id]]$p_val_adj < 0.005]
        pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_logFC > 0 & de.list[[id]]$p_val_adj < 0.005]
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
             go.up.list[help_sort_func(id)] = list(NULL)
             print("null")
             next
        }

        go <- enrichGO(genes_e$ENTREZID, OrgDb = orgdb, ont='ALL',pAdjustMethod = 'BH',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID', readable=TRUE)
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
        genes = de.list[[id]]$gene[de.list[[id]]$avg_logFC < 0 & de.list[[id]]$p_val_adj < 0.005]
        pvals = de.list[[id]]$p_val_adj[de.list[[id]]$avg_logFC < 0 & de.list[[id]]$p_val_adj < 0.005]
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
             go.down.list[help_sort_func(id)] = list(NULL)
             print("null")
             next
        }


        go <- enrichGO(genes_e$ENTREZID, OrgDb = orgdb, ont='ALL',pAdjustMethod = 'BH',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID', readable=TRUE)
        #barplot(go,showCategory=20, drop=T, title = paste("GO analysis for LSK Down genes for cluster", id, sep = " "))
        go.down.list[[id]] = go
    }   
    if(length(go.down.list) == 0){
        return(list())
    }
    go.down.list <- go.down.list[sapply(go.down.list, function(x){!is.null(x)})]
    return(go.down.list)
}



if(!interactive()) {
    conf_main()
}

