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


#Analysis_phases
#1. scrna_phase_preprocess
#2. scrna_phase_clustering
#3. scrna_phase_comparing


### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load
conf = c(
       scrna_phase_preprocess     = 0, ## quality check and preprocessing before integration
       scrna_phase_clustering     = 0, ## integration & clustering
       scrna_phase_comparing      = 2, ## DE GO pathway analysis etc. All rest calculating will be stored here
       scrna_cluster_annotation   = 0, ## Annotate clusters according to `cluster_annotation`
       scrna_clusterwise_xcell    = 0, ## remove cells of each cluster according distinct criterion
       scrna_del_mitogenes        = 0, ## !!!DANGEROUS, once deleted, never recovered!!!
       scrna_merge_clusters       = 1, ## merge clusters
       scrna_remove_clusters      = 0, ## remove clusters
       scrna_remove_recluster     = 0) ## remove clusters and recluster with default resolution

### ----------specific settings for some functions ----------------



## name[new cluster name], value[which clusters to merge together]
scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)

