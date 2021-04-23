# scRNA pipeline

### 1. Overview

This pipeline is a standard pipeline for scRNA analysis using Seurat 4. It enable us to run the analysis automatically.  With only tiny settings, the script will do everything for you.

The following is a brief introduction to the file:

**Makefile** show how to run scripts and produce visualization filess.

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
### --------------Initial info----------------------------
PROJECT = "Mouse Blood project" ## set project name
ORGAN = 'Blood'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Bone-Marrow" #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/


# filtering params when create seurat object
MINCELLS  = 5
MINGENES  = 50


INTEGRATION_OPTION = "seurat" ### or harmony

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv"

data_src = c(
      A_MxCre    =   "data/A_MxCre",
      B_MxCre    =   "data/B_MxCre",
      C_Csnk     =   "data/C_Csnk",
      D_Csnk     =   "data/D_Csnk"
)



##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        A_MxCre      =   "MxCre",
        B_MxCre      =   "MxCre",
        C_Csnk       =   "Csnk",
        D_Csnk       =   "Csnk"
)

## Phase_1, set 1 to regressout
preprocess_regressout = c("mito"      = 1,
                          "ribo"      = 0,
                          "cellcycle" = 1)



### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load
conf = c(
       scrna_phase_preprocess     = 2, ## quality check and preprocessing before integration
       scrna_phase_clustering     = 2, ## integration & clustering
       scrna_phase_comparing      = 2, ## DE GO pathway analysis etc. All rest calculating will be stored here
       scrna_cluster_annotation   = 1, ## Annotate clusters according to `cluster_annotation`
       scrna_clusterwise_xcell    = 0, ## remove cells of each cluster according distinct criterion
       scrna_del_mitogenes        = 0, ## !!!DANGEROUS, once deleted, never recovered!!!
       scrna_merge_clusters       = 0, ## merge clusters
       scrna_remove_clusters      = 0, ## remove clusters
       scrna_remove_recluster     = 0) ## remove clusters and recluster with default resolution



### ----------specific settings for some functions ----------------

## name[your operation name], value[dataframe which cluster, percentage to keep]
scrna_clusterwise_filtercell_settings <- list(
  "mito_cluster0,3,4,5"     =  data.frame(type="mito", max_pct=4, min_pct=0, cluster=c(0,3,4,5)),
  "ribo_cluster2,7_filter"  =  data.frame(type="ribo", max_pct=30, min_pct=0, cluster=c(2,7)),
  "mito_cluster6_filter"    =  data.frame(type="mito", max_pct=3, min_pct=0, cluster=6),
  "mito_cluster8_filter"    =  data.frame(type="mito", max_pct=5, min_pct=0, cluster=8)
)


scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)


scrna_remove_clusters = c(1, 3, 6)
scrna_remove_recluster = c(1, 3, 6)

### cluster annotation
from_cluster_slot = "removed_clusters"
cluster_annotation <- c(
    "0" = "Vascular endothelial",
    "1" = "Fibroblasts 1",
    "2" = "Cardiomyocytes 1",
    "3" = "Endothelial Cells 1",
    "4" = "Macrophages",
    "6" = "Pericytes 1",
    "7" = "Cardiomyocytes 2",
    "8" = "Fibroblasts 2",
    "9" = "Cardiomyocytes 3",
    "10" = "Fibroblasts 3",
    "11" = "Lymphatic endothelial",
    "12" = "VSMCs",
    "13" = "Mesothelial cells",
    "14" = "Cardiomyocytes 4",
    "15" = "Lymphocytes",
    "16" = "T-cells 1",
    "18" = "Endothelial cells 2",
    "19" = "Pericytes 2",
    "20" = "T-cells 2",
    "21" = "Pericytes 2",
    "22" = "T-cells 3"
)
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

run_viz_example.sh

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
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 8

Matrix products: default
BLAS/LAPACK: /data/sz753404/miniconda3/envs/Seurat3/lib/libopenblasp-r0.3.12.so

locale:
 [1] LC_CTYPE=zh_CN.utf8        LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] urltools_1.7.3         reshape2_1.4.4         knitr_1.30             crayon_1.4.1           ComplexHeatmap_2.6.0
 [6] EnhancedVolcano_1.8.0  ggrepel_0.9.0          plotly_4.9.3           digest_0.6.27          openxlsx_4.2.3
[11] cowplot_1.1.1          gridExtra_2.3          kableExtra_1.3.1       msigdbr_7.2.1          ReactomePA_1.34.0
[16] org.Hs.eg.db_3.12.0    org.Mm.eg.db_3.12.0    AnnotationDbi_1.52.0   IRanges_2.24.0         S4Vectors_0.28.0
[21] Biobase_2.50.0         BiocGenerics_0.36.0    clusterProfiler_3.18.0 scHCL_0.1.1            scMCA_0.2.0
[26] genesorteR_0.4.3       SoupX_1.5.0            celda_1.6.1            doParallel_1.0.16      iterators_1.0.13
[31] foreach_1.5.1          stringr_1.4.0          data.table_1.13.6      Matrix_1.2-18          clustree_0.4.3
[36] ggraph_2.0.4           WriteXLS_6.0.0         future.apply_1.6.0     future_1.21.0          progeny_1.12.0
[41] Seurat_3.2.3           dplyr_1.0.4            futile.logger_1.4.3    optparse_1.6.6         glue_1.4.2
[46] Hmisc_4.4-1            ggplot2_3.3.3          Formula_1.2-4          survival_3.2-7         lattice_0.20-41

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3              MCMCprecision_0.4.0         scattermore_0.7             tidyr_1.1.2
  [5] bit64_4.0.5                 irlba_2.3.3                 DelayedArray_0.16.0         rpart_4.1-15
  [9] RCurl_1.98-1.2              generics_0.1.0              lambda.r_1.2.4              RSQLite_2.2.1
 [13] shadowtext_0.0.7            RANN_2.6.1                  combinat_0.0-8              bit_4.0.4
 [17] enrichplot_1.10.0           webshot_0.5.2               xml2_1.3.2                  spatstat.data_2.1-0
 [21] httpuv_1.5.5                SummarizedExperiment_1.20.0 assertthat_0.2.1            viridis_0.5.1
 [25] xfun_0.20                   evaluate_0.14               promises_1.1.1              fansi_0.4.2
 [29] assertive.files_0.0-2       igraph_1.2.6                DBI_1.1.0                   htmlwidgets_1.5.3
 [33] purrr_0.3.4                 ellipsis_0.3.1              backports_1.2.1             deldir_0.2-3
 [37] MatrixGenerics_1.2.0        vctrs_0.3.6                 SingleCellExperiment_1.12.0 Cairo_1.5-12.2
 [41] ROCR_1.0-11                 abind_1.4-5                 RcppEigen_0.3.3.9.1         withr_2.4.1
 [45] ggforce_0.3.2               triebeard_0.3.0             checkmate_2.0.0             sctransform_0.3.2
 [49] getopt_1.20.3               mclust_5.4.6                goftest_1.2-2               cluster_2.1.0
 [53] DOSE_3.16.0                 lazyeval_0.2.2              pkgconfig_2.0.3             labeling_0.4.2
 [57] tweenr_1.0.1                GenomeInfoDb_1.26.0         vipor_0.4.5                 nlme_3.1-150
 [61] nnet_7.3-14                 rlang_0.4.10                globals_0.14.0              lifecycle_1.0.0
 [65] miniUI_0.1.1.1              downloader_0.4              extrafontdb_1.0             enrichR_3.0
 [69] rsvd_1.0.3                  ggrastr_0.2.1               polyclip_1.10-0             matrixStats_0.57.0
 [73] lmtest_0.9-38               graph_1.68.0                zoo_1.8-8                   beeswarm_0.2.3
 [77] base64enc_0.1-3             GlobalOptions_0.1.2         ggridges_0.5.2              pheatmap_1.0.12
 [81] png_0.1-7                   viridisLite_0.3.0           rjson_0.2.20                bitops_1.0-6
 [85] KernSmooth_2.23-18          blob_1.2.1                  shape_1.4.5                 qvalue_2.22.0
 [89] parallelly_1.22.0           jpeg_0.1-8.1                gridGraphics_0.5-1          reactome.db_1.74.0
 [93] scales_1.1.1                graphite_1.36.0             memoise_1.1.0               magrittr_2.0.1
 [97] plyr_1.8.6                  ica_1.0-2                   zlibbioc_1.36.0             compiler_4.0.3
[101] scatterpie_0.1.5            ash_1.0-15                  RColorBrewer_1.1-2          clue_0.3-57
[105] fitdistrplus_1.1-3          XVector_0.30.0              listenv_0.8.0               patchwork_1.1.1
[109] pbapply_1.4-3               htmlTable_2.1.0             formatR_1.7                 MASS_7.3-53
[113] mgcv_1.8-33                 tidyselect_1.1.0            MAST_1.16.0                 stringi_1.5.3
[117] proj4_1.0-10                GOSemSim_2.16.1             assertive.numbers_0.0-2     latticeExtra_0.6-29
[121] fastmatch_1.1-0             tools_4.0.3                 circlize_0.4.11             rstudioapi_0.13
[125] foreign_0.8-80              assertive.types_0.0-3       farver_2.1.0                Rtsne_0.15
[129] rvcheck_0.1.8               BiocManager_1.30.10         shiny_1.6.0                 Rcpp_1.0.6
[133] GenomicRanges_1.42.0        ggalt_0.4.0                 later_1.1.0.1               RcppAnnoy_0.0.18
[137] httr_1.4.2                  assertive.properties_0.0-4  colorspace_2.0-0            rvest_0.3.6
[141] tensor_1.5                  reticulate_1.18             splines_4.0.3               uwot_0.1.10
[145] spatstat.utils_2.1-0        graphlayouts_0.7.1          shinythemes_1.1.2           xtable_1.8-4
[149] jsonlite_1.7.2              futile.options_1.0.1        assertive.base_0.0-9        spatstat_1.64-1
[153] tidygraph_1.2.0             R6_2.5.0                    pillar_1.5.1                htmltools_0.5.1.1
[157] mime_0.10                   fastmap_1.1.0               DT_0.17                     BiocParallel_1.24.1
[161] codetools_0.2-18            maps_3.3.0                  fgsea_1.16.0                utf8_1.1.4
[165] tibble_3.1.0                ggbeeswarm_0.6.0            multipanelfigure_2.1.2      leiden_0.3.6
[169] magick_2.5.2                Rttf2pt1_1.3.8              zip_2.1.1                   GO.db_3.12.1
[173] rmarkdown_2.6               munsell_0.5.0               GetoptLong_1.0.4            DO.db_2.9
[177] GenomeInfoDbData_1.2.4      gtable_0.3.0                extrafont_0.17
```
