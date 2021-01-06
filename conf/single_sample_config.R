### --------------Initial info----------------------------
PROJECT = "Mouse Blood project" ## set project name
ORGAN = 'Blood'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Bone-Marrow" #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/


# filtering params when create seurat object
MINCELLS  = 5
MINGENES  = 50


### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv"

data_src = c(
      A_MxCre    =   "data/A_MxCre"
)

##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        A_MxCre      =   "MxCre"
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
  scrna_phase_singleton      = 1
) 

