
# make Rmd report
args <- commandArgs(TRUE)

PROJ = args[1]
CLUSTER = args[2]
SAVE_DIR = args[3]
FUNCS = args[-c(1:3)]

# TODO better param processing
if(grepl("--proj",PROJ,fixed=TRUE)) PROJ=gsub("--proj=","",PROJ,fixed=TRUE)
if(grepl("--cluster",CLUSTER,fixed=TRUE)) CLUSTER=gsub("--cluster=","",CLUSTER,fixed=TRUE)
if(grepl("--save_dir",SAVE_DIR,fixed=TRUE)) SAVE_DIR=gsub("--save_dir=","",SAVE_DIR,fixed=TRUE)


# User needs to have .set_R_libs.R in their home directory
# otherwise .set_R_libs.R in pipeline folder is used
home_dir = file.path("","home",Sys.info()["user"])
if(file.exists(file.path(home_dir,".set_R_libs.R"))){
  source(file.path(home_dir,".set_R_libs.R"))
}else{
  source(".set_R_libs.R")
}

rmarkdown::render(
  'scrna_pipeline_report.Rmd',
  output_file=file.path(paste0('report',PROJ),'scrna_report.html'),
  clean=TRUE,
  params=list(
    cluster=CLUSTER,
    project=PROJ,
    savedir=SAVE_DIR,
    funcs=FUNCS
  )
)
