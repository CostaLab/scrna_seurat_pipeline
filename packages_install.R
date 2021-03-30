

args <- commandArgs(trailingOnly=TRUE)
# TODO add better argparser
if(length(args)==0){
  force_install = FALSE
}else{
  force_install = as.logical(args[1])
  rlibs_folder = as.character(args[2])

  # If rlibs arg is not an, add path to .libPaths
  if(!is.na(rlibs_folder)){
    if(!dir.exists(rlibs_folder)){
      warning("Rlibs dir provided doesn't exist and will be created!")
    }
    dir.create(rlibs_folder, recursive = TRUE)
    .libPaths(c(rlibs_folder, .libPaths()))
    # dir.exists(file.path("","home",Sys.info()["user"]))
  }
}

# System libraries that might be necessary:
# libfontconfig1-dev
# ibharfbuzz-dev
# libfribidi-dev
# libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# libcairo2-dev
# libxt-dev
# gcc-7
# libgfortran4
# libproj-dev

packages_to_install = c(
  "optparse", "futile.logger", "Seurat", "dplyr",
  "future.apply", "WriteXLS", "clustree", "Matrix",
  "data.table", "ggplot2", "Hmisc", "foreach",
  "devtools", "doParallel", "glue", "openxlsx",
  "rmarkdown", "reshape2", "circlize", "BiocManager",
  "kableExtra", "assertthat", "mclust", "shinythemes",
  "systemfonts", "igraph", "proj4", "Cairo", "ggalt",
  "urltools", "downloadthis", "SoupX"
)

bioc_pkgs_to_install = c(
  "org.Mm.eg.db", "clusterProfiler",
  "org.Hs.eg.db", "ComplexHeatmap",
  "EnhancedVolcano","ReactomePA",
  "msigdbr","limma", "celda"
)

github_pkgs_to_install = c(
  "mahmoudibrahim/genesorteR", "ggjlab/scMCA", "immunogenomics/harmony", "ggjlab/scHCL"
)
gith_pkg_names = gsub("[a-z]+/","",github_pkgs_to_install)




# base R
if(!force_install){
  packages_to_install = packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
  bioc_pkgs_to_install = bioc_pkgs_to_install[!(bioc_pkgs_to_install %in% installed.packages()[,"Package"])]
  github_pkgs_to_install = github_pkgs_to_install[!(gith_pkg_names %in% installed.packages()[,"Package"])]
}

if(length(packages_to_install)){
  install.packages(
    pkgs=packages_to_install,
    repos="https://cloud.r-project.org/",
    dependencies = TRUE
  )
}

if(length(bioc_pkgs_to_install)){
  # BioConductor
  library("BiocManager")
  BiocManager::install(
    pkgs=bioc_pkgs_to_install,
    dependencies=TRUE,
    ask=FALSE,
    update=FALSE
  )
}

if(length(github_pkgs_to_install)){
  library("devtools")
  library("withr")

  with_libpaths(
    .libPaths()[1],
    install_github(
      github_pkgs_to_install,
      dependencies=TRUE,
      upgrade=FALSE
    )
  )
}

## install python module rpy2
ret <- system("pip install rpy2")
if (ret != 0) system("pip show rpy2")

## install python module grip(markdown to html)
ret <- system("pip install grip")
if (ret != 0) system("pip show grip")

ret <- system("pip install Jinja2")
if (ret != 0) system("pip show Jinja2")


failed_to_install = packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
failed_to_install_bioc = bioc_pkgs_to_install[!(bioc_pkgs_to_install %in% installed.packages()[,"Package"])]
failed_to_install_gith = gith_pkg_names[!(gith_pkg_names %in% installed.packages()[,"Package"])]

if(length(failed_to_install)){
  warning(paste0("The following R packages failed to install: ", paste0("\"",failed_to_install,"\"",collapse = ", ")))
}

if(length(failed_to_install_bioc)){
  warning(paste0("The following R BioConductor packages failed to install: ", paste0("\"",failed_to_install_bioc,"\"",collapse = ", ")))
}

if(length(failed_to_install_gith)){
  warning(paste0("The following R GitHub packages failed to install: ", paste0("\"",failed_to_install_gith,"\"",collapse = ", ")))
}

if(!length(c(failed_to_install,failed_to_install_bioc,failed_to_install_gith))){
  message("All R packages were installed successfully!")
}
