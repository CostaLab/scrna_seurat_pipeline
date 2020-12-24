R_install_list <- c("optparse", "futile.logger", "Seurat", "dplyr",
                    "future.apply", "WriteXLS", "clustree", "Matrix",
                    "data.table", "ggplot2", "Hmisc", "foreach",
                    "devtools", "doParallel", "glue", "openxlsx",
                    "rmarkdown", "reshape2", "circlize", "BiocManager",
                    "kableExtra", "assertthat", "configr")



Bio_install_list <- c("org.Mm.eg.db", "clusterProfiler",
                      "org.Hs.eg.db", "ComplexHeatmap",
                      "EnhancedVolcano","ReactomePA",
                      "msigdbr", "harmony")

devtools_install_list <- c("mahmoudibrahim/genesorteR", "ggjlab/scMCA")


failed_vec <- c()
sucess_vec <- c()
for(a_package in R_install_list){
    if (!requireNamespace(a_package, quietly = TRUE)){
        install.packages(a_package, repos = "http://cran.us.r-project.org")
        if(!requireNamespace(a_package, quietly = TRUE)){
            failed_vec <- c(failed_vec, a_package)
        }else{
            sucess_vec <- c(sucess_vec, a_package)
        }
    }

}

for(a_package in Bio_install_list){
    if (!requireNamespace(a_package, quietly = TRUE)){
        ret <- BiocManager::install(a_package)
        if(!requireNamespace(a_package, quietly = TRUE)){
            failed_vec <- c(failed_vec, a_package)
        }else{
            sucess_vec <- c(sucess_vec, a_package)
        }
    }
}
for(a_package in devtools_install_list){
    load_name <- strsplit(a_package, "/")[[1]][2]
    if (!requireNamespace(load_name, quietly = TRUE)){
        ret <- devtools::install_github(a_package)
        if(!requireNamespace(load_name, quietly = TRUE)){
            failed_vec <- c(failed_vec, a_package)
        }else{
            sucess_vec <- c(sucess_vec, a_package)
        }
    }
}



## install python module rpy2
ret <- system("pip install rpy2")
if (ret == 0){
    sucess_vec <- c(sucess_vec, "PYTHON:rpy2")
}else{
    failed_vec <- c(failed_vec, "PYTHON:rpy2")
}

## install python module grip(markdown to html)
ret <- system("pip install grip")
if (ret == 0){
    sucess_vec <- c(sucess_vec, "PYTHON:grip")
}else{
    failed_vec <- c(failed_vec, "PYTHON:grip")
}

ret <- system("pip install jinja2")
if (ret == 0){
    sucess_vec <- c(sucess_vec, "PYTHON:jinja2")
}else{
    failed_vec <- c(failed_vec, "PYTHON:jinja2")
}



if(length(failed_vec) == 0){
    message("\n\n===================================================\n
               Congratulations!!!
               Environment have satisfied the pipeline!\n
               \r===================================================\n")
    if(length(sucess_vec) !=0){
        message("succeed in installing packages: ", paste(sucess_vec, collapse=" "))
    }
}else{
    message("\n\n===================================================\n
             Failed to install packages: ", paste(failed_vec, collapse=" "), " please check!!!!
            \r===================================================\n")
}
