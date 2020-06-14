R_install_list <- c("optparse", "futile.logger", "Seurat", "dplyr", 
                    "future.apply", "WriteXLS", "clustree", "Matrix", 
                    "data.table", "ggplot2", "Hmisc", "foreach", "doParallel")



Bio_install_list <- c("scMCA", "clusterProfiler", "org.Mm.eg.db", 
                      "clusterProfiler", "org.Mm.eg.db")


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


if(length(failed_vec) == 0){
    message("\n\n===================================================\n
               Congratulations!!!
               Environment have satisfied the pipeline!\n
               \r===================================================\n")
    if(length(sucess_vec) !=0){
        message("succeed in installing packages:", paste(sucess_vec))    
    } 
}else{
    message("\n\n===================================================\n
              Failed to install packages", paste(failed_vec), " please check!!!!
            \r===================================================\n")
}
