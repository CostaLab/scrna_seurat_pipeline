
suppressPackageStartupMessages(library(archive))

save_object <- function(object, file_name, file_format=NULL){

  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))

  if(file_format %in% "nocomp"){
    saveRDS(object = object, file = file_name, compress = FALSE)
    return(invisible(NULL))
  }

  if(file_format %in% c("zstd", "lz4")){
    con <- archive::file_write(file = file_name, filter = file_format)
    open(con)
    saveRDS(object = object, file = con)
    close(con)
  }else{
    saveRDS(object = object, file = file_name, compress = file_format)
  }
}

load_object <- function(file_name){
  con <- archive::file_read(file = file_name)
  res <- readRDS(file = con)
  close(con)
  return(res)
}


# use_tools: if True, use the tools directory to load the object
#            if False, use the save_dir
seutools_partition <- function(scrna, partition, save_dir, allinone=FALSE, use_tools=FALSE){
  #scrna@tools, store into partition if allinone is FALSE
  out <- NULL
  if(allinone == FALSE){
    #assertthat::assert_that(scrna@tools$allinone == FALSE)
    if(!(partition %in% names(scrna@tools))){
      stop("partition not found in scrna@tools")
    }
    if(endsWith(scrna@tools[[partition]], "Rds")){
      if(use_tools == TRUE){
        out <- load_object(scrna@tools[[partition]])
      }else{
        out <- load_object(file.path(save_dir, "partition", glue::glue("{partition}.Rds")))
      }
    }else{
      stop(glue::glue("The {partition}.Rds is not existing!"))
    }
  }else{
    out <- scrna@tools[[partition]]
  }

  return(out)
}



seu_assay <- function(scrna, assay, save_dir, allinone=FALSE, use_tools=FALSE){
  #scrna@tools, store into assay if allinone is FALSE
  out <- NULL
  if(allinone == FALSE){
    #assertthat::assert_that(scrna@tools$allinone == FALSE)
    if(!((assay %in% names(scrna@tools$assay_info)) | (assay %in% names(scrna@assays)))){
      stop("assay not found in scrna@tools$assay_info or scrna@assay")
    }
    if(endsWith(scrna@tools$assay_info[[assay]], "Rds")){
      if(use_tools == TRUE){
        assay_data <- load_object(scrna@tools$assay_info[[assay]])
        assertthat::assert_that(all(colnames(scrna) %in% colnames(assay_data$assay)))
        scrna[[assay]] <- subset(assay_data$assay, cells=colnames(scrna))
        rm(assay_data)
      }else{
        assay_data <- load_object(file.path(save_dir, "assays", glue::glue("{assay}.Rds")))
        assertthat::assert_that(all(colnames(scrna) %in% colnames(assay_data$assay)))
        scrna[[assay]] <- subset(assay_data$assay, cells=colnames(scrna))
        rm(assay_data)
      }
    }else{
      stop(glue::glue("The {assay}.Rds is not existing!"))
    }
  }
  return(scrna)
}
