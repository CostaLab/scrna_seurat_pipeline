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
|   |-- config_remove_reclusters.R
|   |-- config_removed_clusters.R
|   `-- single_sample_config.R
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
    |-- 3_external_markers.Rmd
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

     NK1_Gli1_IRI     =   "data/summed_mtx/sum_NK1_Gli1_IRI/",
     NK2_CD45_IRI     =   "data/summed_mtx/sum_NK2_CD45_IRI/",
     NK3_Gli1_Sham    =   "data/summed_mtx/sum_NK3_Gli1_Sham/",
     NK4_CD45_Sham    =   "data/summed_mtx/sum_NK4_CD45_Sham/"
)

#A_MxCre B_MxCre  C_Csnk  D_Csnk
##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
    NK1_Gli1_IRI      = "IRI",
    NK2_CD45_IRI      = "IRI",
    NK3_Gli1_Sham     = "Sham",
    NK4_CD45_Sham     = "Sham"
)


### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load   3.
conf = c(
       scrna_rawdata              = 1, ## read count matrix and merge samples to a Seurat obj
       scrna_filter               = 1, ## filter nFeatureRNA and nCountRNA
       scrna_preprocess           = 1, ## Normalize & FindvariablegFeatures and ScaleData
       scrna_cellcycle            = 1, ## Cell cycle scoring
       scrna_cycleRegressOut      = 1, ## Regress out cell cycle effects
       scrna_CCmitoRegressOut     = 1, ## Regress out mito & cell cycle
       scrna_mitoRegressOut       = 0, ## Regress out mito only
       scrna_riboRegressOut       = 0, ## Regress out ribo only
       scrna_CCriboRegressOut     = 0, ## Regress out ribo & cell cycle
       scrna_mitoRiboRegressOut   = 0, ## Regress out ribo&mito
       scrna_CCmitoRiboRegressOut = 0, ## Regress out ribo&mito  & cell cycle
       scrna_integration          = 1, ## Integrate samples using Seurat 3
       scrna_ScaleIntegration     = 1, ## ScaleDatai&PCA and UMAP
       scrna_batchclustering      = 1, ## clustering with resolution from 0.1 to 0.8
       scrna_batch_markergenes    = 1, ## Marker Genes for clusters with different resolutions
       scrna_clustering           = 1, ## Set seurat_clusters or re-calculate
       scrna_clusterwise_xcell    = 1, ## remove cells of each cluster according distinct criterion
       scrna_del_mitogenes        = 0, ## !!!DANGEROUS, once deleted, never recovered!!!
       scrna_markergenes          = 1, ## markergenes for seurat_clusters
       scrna_genesorteR           = 1, ## genesorteR analysis
       scrna_go                   = 1, ## Gene Ontology analysis
       scrna_kegg                 = 1, ## kegg enrichment analysis
       scrna_reactome             = 1, ## reactome enrichment analysis
       scrna_hallmark             = 1, ## hallmark enrichment analysis
       scrna_fishertest_clusters  = 1, ## fisher test for clusters and stages
       scrna_MCAannotate          = 1, ## scMCA annotation celltypes
       scrna_ExternalAnnotation   = 1, ## Annotation from given databases(tsv)
       scrna_dego_name            = 1, ## DE & GO between samples
       scrna_dego_stage           = 1, ## DE & GO between stages
       scrna_dego_stage_vsRest    = 1, ## DE & GO between one stage and all Rest
       scrna_pathway_name         = 1, ## samples comparison KEGG&Reactome&hallmark
       scrna_pathway_stage        = 1, ## stages comparison KEGG&Reactome&hallmark
       scrna_pathway_stage_vsRest = 1, ## stages vsRest comparison KEGG&Reactome&hallmark
       scrna_clusterwise_xcell    = 0, ## keep cells for each cluster according to mito&ribo
       scrna_fishertest_clusters  = 0, ## fisher test for clusters and stages
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

Sometimes, we have to tune some parameters or to parallelize the program, this is where parameters of **data_factory.R** comes from.  For example,  set *-n* to decide how many cores do you need to run.

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
r$> sessionInfo()
R version 4.0.1 (2020-06-06)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 8 (Core)

Matrix products: default
BLAS/LAPACK: /mnt/data/mingbo/miniconda3/envs/r4.0.1/lib/libopenblasp-r0.3.9.so

locale:
 [1] LC_CTYPE=en_US.UTF-8    LC_NUMERIC=C            LC_TIME=en_US.UTF-8     LC_COLLATE=en_US.UTF-8  LC_MONETARY=C
 [6] LC_MESSAGES=en_US.UTF-8 LC_PAPER=C              LC_NAME=C               LC_ADDRESS=C            LC_TELEPHONE=C
[11] LC_MEASUREMENT=C        LC_IDENTIFICATION=C

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] EnhancedVolcano_1.7.16 ggrepel_0.8.2          ComplexHeatmap_2.6.0   openxlsx_4.2.2         digest_0.6.27          gridExtra_2.3
 [7] urltools_1.7.3         reshape2_1.4.4         cowplot_1.1.0          kableExtra_1.2.1       knitr_1.30             msigdbr_7.2.1
[13] ReactomePA_1.33.1      org.Hs.eg.db_3.12.0    org.Mm.eg.db_3.12.0    AnnotationDbi_1.52.0   IRanges_2.24.0         S4Vectors_0.28.0
[19] Biobase_2.50.0         BiocGenerics_0.36.0    clusterProfiler_3.17.5 scMCA_0.2.0            genesorteR_0.4.3       doParallel_1.0.16
[25] iterators_1.0.13       foreach_1.5.1          Hmisc_4.4-1            Formula_1.2-4          survival_3.2-7         lattice_0.20-41
[31] stringr_1.4.0          data.table_1.13.2      Matrix_1.2-18          clustree_0.4.3         ggraph_2.0.3           ggplot2_3.3.2
[37] WriteXLS_6.0.0         future.apply_1.6.0     future_1.19.1          glue_1.4.2             dplyr_1.0.2            Seurat_3.2.2
[43] futile.logger_1.4.3    optparse_1.6.6

loaded via a namespace (and not attached):
  [1] reticulate_1.16       tidyselect_1.1.0      RSQLite_2.2.1         htmlwidgets_1.5.2     BiocParallel_1.24.0   Rtsne_0.15
  [7] scatterpie_0.1.5      munsell_0.5.0         codetools_0.2-17      ica_1.0-2             DT_0.16               miniUI_0.1.1.1
 [13] withr_2.3.0           colorspace_1.4-1      GOSemSim_2.15.2       ggalt_0.4.0           rstudioapi_0.11       ROCR_1.0-11
 [19] tensor_1.5            Rttf2pt1_1.3.8        DOSE_3.15.4           listenv_0.8.0         polyclip_1.10-0       bit64_4.0.5
 [25] farver_2.0.3          pheatmap_1.0.12       downloader_0.4        vctrs_0.3.4           generics_0.1.0        lambda.r_1.2.4
 [31] xfun_0.18             R6_2.5.0              ggbeeswarm_0.6.0      clue_0.3-57           graphlayouts_0.7.0    rsvd_1.0.3
 [37] spatstat.utils_1.17-0 fgsea_1.15.2          promises_1.1.1        scales_1.1.1          nnet_7.3-14           enrichplot_1.9.4
 [43] beeswarm_0.2.3        gtable_0.3.0          ash_1.0-15            Cairo_1.5-12.2        globals_0.13.1        goftest_1.2-2
 [49] tidygraph_1.2.0       rlang_0.4.8           GlobalOptions_0.1.2   splines_4.0.1         extrafontdb_1.0       lazyeval_0.2.2
 [55] checkmate_2.0.0       BiocManager_1.30.10   abind_1.4-5           backports_1.2.0       httpuv_1.5.4          qvalue_2.21.0
 [61] extrafont_0.17        tools_4.0.1           ellipsis_0.3.1        RColorBrewer_1.1-2    ggridges_0.5.2        Rcpp_1.0.5
 [67] plyr_1.8.6            base64enc_0.1-3       purrr_0.3.4           rpart_4.1-15          deldir_0.1-29         GetoptLong_1.0.4
 [73] pbapply_1.4-3         viridis_0.5.1         zoo_1.8-8             cluster_2.1.0         magrittr_1.5          futile.options_1.0.1
 [79] DO.db_2.9             circlize_0.4.11       triebeard_0.3.0       lmtest_0.9-38         RANN_2.6.1            reactome.db_1.70.0
 [85] fitdistrplus_1.1-1    matrixStats_0.57.0    patchwork_1.0.1       mime_0.9              evaluate_0.14         xtable_1.8-4
 [91] jpeg_0.1-8.1          mclust_5.4.6          shape_1.4.5           compiler_4.0.1        maps_3.3.0            tibble_3.0.4
 [97] KernSmooth_2.23-17    crayon_1.3.4          shadowtext_0.0.7      htmltools_0.5.0       mgcv_1.8-33           later_1.1.0.1
[103] tidyr_1.1.2           DBI_1.1.0             tweenr_1.0.1          formatR_1.7           proj4_1.0-10          MASS_7.3-53
[109] rappdirs_0.3.1        getopt_1.20.3         igraph_1.2.6          pkgconfig_2.0.3       rvcheck_0.1.8         foreign_0.8-80
[115] plotly_4.9.2.1        xml2_1.3.2            vipor_0.4.5           webshot_0.5.2         rvest_0.3.6           sctransform_0.3.1
[121] RcppAnnoy_0.0.16      graph_1.67.1          spatstat.data_1.4-3   rmarkdown_2.4         leiden_0.3.3          fastmatch_1.1-0
[127] htmlTable_2.1.0       uwot_0.1.8            shiny_1.5.0           graphite_1.35.3       rjson_0.2.20          lifecycle_0.2.0
[133] nlme_3.1-149          jsonlite_1.7.1        viridisLite_0.3.0     pillar_1.4.6          ggrastr_0.2.1         fastmap_1.0.1
[139] httr_1.4.2            GO.db_3.12.1          zip_2.1.1             spatstat_1.64-1       png_0.1-7             shinythemes_1.1.2
[145] bit_4.0.4             ggforce_0.3.2         stringi_1.5.3         blob_1.2.1            latticeExtra_0.6-29   memoise_1.1.0
[151] irlba_2.3.3
```
