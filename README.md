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
