
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

