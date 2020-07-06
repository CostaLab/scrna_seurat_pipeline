# scRNA pipeline 

### 1. Overview

This pipeline is a standard pipeline for scRNA analysis using Seurat 3. It enable us to run the analysis automatically.  With only tiny settings, the script will do everything for you. 

The following is a brief introduction to the file:

Two **run.sh** files show how to run scripts and produce visualization filess.

The main pipeline script is **data_factory.R**. Please set conf/config.R before run it. 

The **config.R** file settings include the Count Matrix path, stages map settings and functions you plan to run and a few initial settings. You can also specify your own config.R file according to your needs. 

**packages_install.R** is a simple script allows you to install all packages needed.

Files in **viz** directory is to visualize using the data produced by **data_factory.R**. Some of them are generated from template files

Meanwhile, **code_generator.py** can generate some R markdown files and index.md according to **the config.R**. 

```shell

|-- conf
|   |-- config.R
|   |-- config_merged_clusters.R
|   |-- config_removed_clusters.R
|   `-- config_remove_recluster.R
|-- data_factory.R
|-- external  
|   `-- Human_and_mouse_cell_markers-Markers.tsv
|-- packages_install.R
|-- run.sh
`-- viz
    |-- 1_quality_report.Rmd
    |-- 2_clustering.Rmd
    |-- 2_clusters_DEs.Rmd
    |-- 3_DE_GO-analysis.Rmd
    |-- code_generator.py
    |-- run.sh
    `-- template
        |-- DE-GO-1v1.template
        |-- DE-GO-stagesVS.template
        `-- index.template
```



### 2. How it works

#### conf/config.R

The pipeline is controlled by the most important configuration file **conf/config.R**. To configure *data_src* to tell the script where are the Count Matrix located.  And *stage_list* set several groups that may be useful in the downstream analysis. *conf*  has defined a bunch of functions which have been implemented in the main script **data_factory.R **.  The setting is really flexible for customizing. Set to 2 to load the existed data generated from the function while set to 1 to run the function. Set it to 0 if you're not going to run the function. For example, for some datasets, there's no need to remove the cell cycle effect,  just set it to 0.  The following is an example settings.

```R
#conf/config.R
### --------------Initail info----------------------------
PROJECT = "Intestine project" ## set project name 
ORGAN = 'Intestine'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Fetal_Intestine"  #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/ 

MINCELLS  = 5   ## when creating seuratobject
MINGENES  = 50  ## when creating seuratobject

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv" 

data_src = c( 
      I       =    "./results.200402/I/filtered_feature_bc_matrix",
      Nd7     =    "./results.200402/Nd7/filtered_feature_bc_matrix",
      Ad7     =    "./results.200402/Ad7/filtered_feature_bc_matrix"
)

#A_MxCre B_MxCre  C_Csnk  D_Csnk 
##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        I      =   "I",
        Nd7    =   "Nd7",
        Ad7    =   "Ad7"
)

### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load   3.  
conf = c(
       scrna_rawdata              = 1, ## read count matrix and merge samples to a Seurat obj 
       scrna_filter               = 1, ## filter nFeatureRNA and nCountRNA
       scrna_preprocess           = 1, ## Normalize & FindvariablegFeatures and ScaleData
       scrna_excludeRegressOut    = 1, ## Regress out mito and ribo effects 
       scrna_cellcycle            = 1, ## Cell cycle scoring
       scrna_cycleRegressOut      = 1, ## Regress out cell cycle effects
       scrna_RegressOutAll        = 1, ## Regress out celly cycle & mito & ribo
       scrna_integration          = 1, ## Integrate samples using Seurat 3
       scrna_ScaleIntegration     = 1, ## ScaleDatai&PCA and UMAP
       scrna_batchclustering      = 1, ## clustering with resolution from 0.1 to 0.8
       scrna_batch_markergenes    = 1, ## Marker Genes for clusters with different resolutions
       scrna_clustering           = 1, ## Set seurat_clusters or re-calculate
       scrna_markergenes          = 1, ## markergenes for seurat_clusters
       scrna_go                   = 1, ## Gene Ontology analysis
       scrna_fishertest_clusters  = 1, ## fisher test for clusters and stages
       scrna_MCAannotate          = 1, ## scMCA annotation celltypes
       scrna_ExternalAnnotation   = 1, ## Annotation from given databases(tsv)
       scrna_dego_name            = 1, ## DE & GO between samples
       scrna_dego_stage           = 1, ## GO down for mark genes
       scrna_merge_clusters       = 0, ## merge clusters 
       scrna_remove_clusters      = 0, ## remove clusters 
       scrna_remove_recluster     = 0, ## remove clusters and recluster with defualt resolution
       scrna_markergenes          = 0, ## markergenes for seurat_clusters
       scrna_go                   = 0, ## Gene Ontology analysis
       scrna_dego_name            = 0, ## DE & GO between samples
       scrna_dego_stage           = 0) ## GO down for mark genes




### ----------specific settings for some functions ----------------

scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)


scrna_remove_clusters = c(1, 3, 6)
scrna_remove_recluster = c(1, 3, 6)
```

#### Parameters

Sometimes, we have to tune some parameters or to parallelize the program, this is where parameters of **data_factory.R ** comes from.  For example,  set *-n* to decide how many cores do you need to run.  

```shell
Rscript data_factory.R --help
```



### 3. Examples
In this section, we show several examples how to run the programme with some paramters. For example, to do the analysis with removed clusters. We have placed configuration examples to conf directories that you can use(removed_cluster, etc.). Notice that you have to specify which clusters to use in this case. By default, the option is seurat_clusters. For other the downstream analysis use removed_clusters, remove_recluster, merged_clusters.



#### Produce Rdata

```shell
## 50 cores run, future parallel memory
Rscript data_factory.R -n 50  --MaxMemMega=100000

# set config file to use, default(conf/config.R) 
Rscript data_factory.R -c conf/another_config.R 


# set count matrix format, 10X and 10X_h5 
Rscript data_factory.R -f 10X_h5 


# set cluster resolution(0.5) for metadata slot "seurat_clusters"
Rscript data_factory.R -r 0.5


# set Find anchors(usually cca) for integration(1:20), dims to use
Rscript data_factory.R -d 1:20 


# set Find neighbours for(usaully pca after integration) clustering(1:12), dims to use
Rscript data_factory.R -x 1:12


# set using clusters: "merged_clusters". Do DE or GO analysis
Rscript data_factory.R -a merged_clusters -c conf/config_merged_clusters.R


# set using clusters: "removed_clusters". Do DE or GO analysis
Rscript data_factory.R -a removed_clusters -c conf/config_remove_clusters.R


# set using clusters: "remove_recluster". Do DE or GO analysis
Rscript data_factory.R -a remove_recluster -c conf/config_remove_recluster.R
```




#### visualisation

```shell
## set a parameter when visualisation 

viz/run.sh

# In general, we choose seurat_clusters, 
# If you are using removed or merged clusters,
# choose the following:
            # seurat_clusters
            # merged_clusters
            # removed_clusters
            # remove_recluster 

## Here we set normal senario clusters
cluster="seurat_clusters" # <seurat_clusters|removed_clusters|remove_recluster|merged_recluster>
```



#### R sessionInfo 

```R
> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 8 (Core)

Matrix products: default
BLAS/LAPACK: /mnt/data/mingbo/anaconda3/envs/rs3/lib/libopenblasp-r0.3.9.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] org.Mm.eg.db_3.10.0    AnnotationDbi_1.48.0   IRanges_2.20.2
 [4] S4Vectors_0.24.4       Biobase_2.46.0         BiocGenerics_0.32.0
 [7] clusterProfiler_3.14.3 scMCA_0.2.0            doParallel_1.0.15
[10] iterators_1.0.12       foreach_1.5.0          Hmisc_4.4-0
[13] Formula_1.2-3          survival_3.1-11        lattice_0.20-41
[16] stringr_1.4.0          data.table_1.12.8      Matrix_1.2-18
[19] clustree_0.4.2         ggraph_2.0.2           ggplot2_3.3.0
[22] WriteXLS_5.0.0         future.apply_1.4.0     future_1.16.0
[25] dplyr_0.8.5            Seurat_3.1.4           futile.logger_1.4.3
[28] optparse_1.6.4

loaded via a namespace (and not attached):
  [1] reticulate_1.15      tidyselect_1.0.0     RSQLite_2.2.0
  [4] htmlwidgets_1.5.1    grid_3.6.2           BiocParallel_1.20.1
  [7] Rtsne_0.15           munsell_0.5.0        codetools_0.2-16
 [10] mutoss_0.1-12        ica_1.0-2            DT_0.13
 [13] withr_2.2.0          colorspace_1.4-1     GOSemSim_2.12.1
 [16] knitr_1.28           rstudioapi_0.11      ROCR_1.0-7
 [19] DOSE_3.12.0          gbRd_0.4-11          listenv_0.8.0
 [22] Rdpack_0.11-1        urltools_1.7.3       mnormt_1.5-6
 [25] polyclip_1.10-0      bit64_0.9-7          farver_2.0.3
 [28] pheatmap_1.0.12      vctrs_0.2.4          TH.data_1.0-10
 [31] lambda.r_1.2.4       xfun_0.13            R6_2.4.1
 [34] graphlayouts_0.6.0   rsvd_1.0.3           gridGraphics_0.5-0
 [37] fgsea_1.12.0         bitops_1.0-6         assertthat_0.2.1
 [40] promises_1.1.0       scales_1.1.0         enrichplot_1.6.1
 [43] multcomp_1.4-12      nnet_7.3-13          gtable_0.3.0
 [46] npsurv_0.4-0         globals_0.12.5       tidygraph_1.1.2
 [49] sandwich_2.5-1       rlang_0.4.6          splines_3.6.2
 [52] lazyeval_0.2.2       acepack_1.4.1        europepmc_0.3
 [55] checkmate_2.0.0      BiocManager_1.30.10  reshape2_1.4.3
 [58] backports_1.1.7      httpuv_1.5.2         qvalue_2.18.0
 [61] tools_3.6.2          ggplotify_0.0.5      ellipsis_0.3.0
 [64] gplots_3.0.3         RColorBrewer_1.1-2   ggridges_0.5.2
 [67] TFisher_0.2.0        Rcpp_1.0.4           plyr_1.8.6
 [70] progress_1.2.2       base64enc_0.1-3      purrr_0.3.3
 [73] prettyunits_1.1.1    rpart_4.1-15         pbapply_1.4-2
 [76] viridis_0.5.1        cowplot_1.0.0        zoo_1.8-7
 [79] ggrepel_0.8.2        cluster_2.1.0        magrittr_1.5
 [82] futile.options_1.0.1 DO.db_2.9            triebeard_0.3.0
 [85] lmtest_0.9-37        RANN_2.6.1           mvtnorm_1.1-0
 [88] fitdistrplus_1.0-14  hms_0.5.3            patchwork_1.0.0
 [91] lsei_1.2-0           mime_0.9             xtable_1.8-4
 [94] jpeg_0.1-8.1         gridExtra_2.3        compiler_3.6.2
 [97] tibble_3.0.0         KernSmooth_2.23-16   crayon_1.3.4
[100] htmltools_0.4.0      later_1.0.0          tidyr_1.0.2
[103] DBI_1.1.0            tweenr_1.0.1         formatR_1.7
[106] MASS_7.3-51.5        rappdirs_0.3.1       getopt_1.20.3
[109] cli_2.0.2            gdata_2.18.0         metap_1.3
[112] igraph_1.2.5         pkgconfig_2.0.3      sn_1.6-1
[115] rvcheck_0.1.8        numDeriv_2016.8-1.1  foreign_0.8-76
[118] plotly_4.9.2         xml2_1.3.2           multtest_2.42.0
[121] bibtex_0.4.2.2       digest_0.6.25        sctransform_0.2.1
[124] RcppAnnoy_0.0.16     tsne_0.1-3           fastmatch_1.1-0
[127] leiden_0.3.3         htmlTable_1.13.3     uwot_0.1.8
[130] shiny_1.4.0.2        gtools_3.8.2         lifecycle_0.2.0
[133] nlme_3.1-145         jsonlite_1.6.1       viridisLite_0.3.0
[136] fansi_0.4.1          pillar_1.4.3         GO.db_3.10.0
[139] fastmap_1.0.1        httr_1.4.1           plotrix_3.7-7
[142] glue_1.4.1           png_0.1-7            shinythemes_1.1.2
[145] bit_1.1-15.2         ggforce_0.3.1        stringi_1.4.6
[148] blob_1.2.1           latticeExtra_0.6-29  caTools_1.18.0
[151] memoise_1.1.0        irlba_2.3.3          ape_5.3
```
