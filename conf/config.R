### --------------Initail info----------------------------
PROJECT = "Intestine project" ## set project name 
ORGAN = 'Intestine'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Fetal_Intestine"  #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/ 

# filtering params when create seurat object
MINCELLS  = 5 
MINGENES  = 50

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv" 

data_src = c( 
      I = "./results.200402/I/filtered_feature_bc_matrix",
      Nd7 = "./results.200402/Nd7/filtered_feature_bc_matrix",
      Ad7 = "./results.200402/Ad7/filtered_feature_bc_matrix"
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
