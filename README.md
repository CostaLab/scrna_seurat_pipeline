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

## If genesets you need are not included, please attach your geneset to the gmt.gz file.
MSigDB_GENESET_HUMAN_GMT_FILE  = "external/Human_msigdb.v7.2.symbols.gmt.gz"

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


## Genesets candidate names please check external/MSigDB_names.txt
MSigDB_Geneset_names <- c(
    "NABA_COLLAGENS",
    "NABA_SECRETED_FACTORS",
    "NABA_ECM_GLYCOPROTEINS",
    "NABA_CORE_MATRISOME",
    "NABA_ECM_REGULATORS",
    "NABA_MATRISOME_ASSOCIATED",
    "NABA_ECM_AFFILIATED",
    "NABA_BASEMENT_MEMBRANES",
    "NABA_PROTEOGLYCANS",
    "NABA_MATRISOME"
)


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
R version 4.2.2 (2022-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /data/sz753404/miniconda3/envs/R422/lib/libopenblasp-r0.3.21.so

locale:
[1] C

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] genesorteR_0.4.3            SoupX_1.6.2                 celda_1.14.2                SingleCellExperiment_1.20.0
 [5] SummarizedExperiment_1.28.0 GenomicRanges_1.50.0        GenomeInfoDb_1.34.8         MatrixGenerics_1.10.0      
 [9] matrixStats_0.63.0          doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2              
[13] data.table_1.14.8           Matrix_1.5-3                future.apply_1.10.0         future_1.32.0              
[17] progeny_1.20.0              futile.logger_1.4.3         DOSE_3.24.2                 plotly_4.10.1              
[21] digest_0.6.31               WriteXLS_6.4.0              gridExtra_2.3               ggsci_3.0.0                
[25] ggridges_0.5.4              SeuratObject_4.1.3          Seurat_4.3.0                EnhancedVolcano_1.16.0     
[29] ggrepel_0.9.3               ComplexHeatmap_2.14.0       openxlsx_4.2.5.2            stringr_1.5.0              
[33] clustree_0.5.0              ggraph_2.1.0                urltools_1.7.3              reshape2_1.4.4             
[37] cowplot_1.1.1               dplyr_1.1.0                 kableExtra_1.3.4            knitr_1.42                 
[41] Hmisc_5.0-1                 glue_1.6.2                  crayon_1.5.2                optparse_1.7.3             
[45] msigdbr_7.5.1               ReactomePA_1.42.0           org.Hs.eg.db_3.16.0         org.Mm.eg.db_3.16.0        
[49] clusterProfiler_4.7.1.003   GSEABase_1.60.0             graph_1.76.0                annotate_1.76.0            
[53] XML_3.99-0.13               AnnotationDbi_1.60.0        IRanges_2.32.0              S4Vectors_0.36.0           
[57] Biobase_2.58.0              BiocGenerics_0.44.0         scHCL_0.1.1                 scMCA_0.2.0                
[61] ggplot2_3.4.1              

loaded via a namespace (and not attached):
  [1] ica_1.0-3                  svglite_2.1.1              assertive.properties_0.0-5 lmtest_0.9-40             
  [5] MASS_7.3-58.3              nlme_3.1-162               backports_1.4.1            GOSemSim_2.24.0           
  [9] rlang_1.1.0                XVector_0.38.0             HDO.db_0.99.1              ROCR_1.0-11               
 [13] irlba_2.3.5.1              BiocParallel_1.32.5        rjson_0.2.21               bit64_4.0.5               
 [17] pheatmap_1.0.12            sctransform_0.3.5          spatstat.sparse_3.0-1      spatstat.geom_3.1-0       
 [21] tidyselect_1.2.0           fitdistrplus_1.1-8         tidyr_1.3.0                assertive.types_0.0-3     
 [25] zoo_1.8-11                 xtable_1.8-4               magrittr_2.0.3             evaluate_0.20             
 [29] cli_3.6.0                  zlibbioc_1.44.0            rstudioapi_0.14            miniUI_0.1.1.1            
 [33] sp_1.6-0                   rpart_4.1.19               fastmatch_1.1-3            lambda.r_1.2.4            
 [37] treeio_1.22.0              RcppEigen_0.3.3.9.3        shiny_1.7.4                xfun_0.37                 
 [41] clue_0.3-64                gson_0.1.0                 cluster_2.1.4              tidygraph_1.2.3           
 [45] KEGGREST_1.38.0            tibble_3.2.0               ape_5.7-1                  listenv_0.9.0             
 [49] Biostrings_2.66.0          png_0.1-8                  withr_2.5.0                bitops_1.0-7              
 [53] ggforce_0.4.1              plyr_1.8.8                 assertive.base_0.0-9       pillar_1.8.1              
 [57] GlobalOptions_0.1.2        cachem_1.0.7               GetoptLong_1.0.5           graphite_1.44.0           
 [61] vctrs_0.5.2                ellipsis_0.3.2             generics_0.1.3             tools_4.2.2               
 [65] foreign_0.8-84             munsell_0.5.0              tweenr_2.0.2               fgsea_1.24.0              
 [69] DelayedArray_0.24.0        fastmap_1.1.1              compiler_4.2.2             abind_1.4-5               
 [73] httpuv_1.6.9               GenomeInfoDbData_1.2.9     enrichR_3.1                lattice_0.20-45           
 [77] deldir_1.0-6               utf8_1.2.3                 later_1.3.0                jsonlite_1.8.4            
 [81] multipanelfigure_2.1.2     scales_1.2.1               tidytree_0.4.2             pbapply_1.7-0             
 [85] lazyeval_0.2.2             promises_1.2.0.1           goftest_1.2-3              spatstat.utils_3.0-2      
 [89] reticulate_1.25            checkmate_2.1.0            rmarkdown_2.20             webshot_0.5.4             
 [93] Rtsne_0.16                 downloader_0.4             uwot_0.1.14                igraph_1.4.1              
 [97] survival_3.5-5             systemfonts_1.0.4          htmltools_0.5.4            memoise_2.0.1             
[101] graphlayouts_0.8.4         viridisLite_0.4.1          mime_0.12                  rappdirs_0.3.3            
[105] futile.options_1.0.1       RSQLite_2.3.0              yulab.utils_0.0.6          blob_1.2.3                
[109] shinythemes_1.2.0          splines_4.2.2              Formula_1.2-5              RCurl_1.98-1.10           
[113] assertive.numbers_0.0-2    colorspace_2.1-0           base64enc_0.1-3            shape_1.4.6               
[117] assertive.files_0.0-2      aplot_0.1.10               nnet_7.3-18                mclust_6.0.0              
[121] Rcpp_1.0.10                RANN_2.6.1                 circlize_0.4.15            enrichplot_1.18.3         
[125] fansi_1.0.4                parallelly_1.34.0          R6_2.5.1                   lifecycle_1.0.3           
[129] formatR_1.14               zip_2.2.2                  curl_4.3.3                 leiden_0.4.3              
[133] getopt_1.20.3              qvalue_2.30.0              RcppAnnoy_0.0.20           RColorBrewer_1.1-3        
[137] spatstat.explore_3.1-0     htmlwidgets_1.6.1          polyclip_1.10-4            triebeard_0.4.1           
[141] purrr_1.0.1                shadowtext_0.1.2           gridGraphics_0.5-1         reactome.db_1.82.0        
[145] rvest_1.0.3                globals_0.16.2             htmlTable_2.4.1            patchwork_1.1.2           
[149] spatstat.random_3.1-4      progressr_0.13.0           codetools_0.2-19           GO.db_3.16.0              
[153] MCMCprecision_0.4.0        gtable_0.3.1               DBI_1.1.3                  ggfun_0.0.9               
[157] tensor_1.5                 httr_1.4.5                 KernSmooth_2.23-20         stringi_1.7.12            
[161] farver_2.1.1               viridis_0.6.2              magick_2.7.4               ggtree_3.6.2              
[165] DT_0.27                    xml2_1.3.3                 combinat_0.0-8             ggplotify_0.1.0           
[169] scattermore_0.8            bit_4.0.5                  scatterpie_0.1.8           spatstat.data_3.0-1       
[173] pkgconfig_2.0.3            babelgene_22.9            

```
