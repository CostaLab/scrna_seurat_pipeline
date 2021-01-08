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
       scrna_phase_preprocess     = 1, ## quality check and preprocessing before integration
       scrna_phase_clustering     = 1, ## integration & clustering
       scrna_phase_comparing      = 1, ## DE GO pathway analysis etc. All rest calculating will be stored here
       scrna_cluster_annotation   = 0, ## Annotate clusters according to `cluster_annotation`
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

## name[new cluster name], value[which clusters to merge together]
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



### MCA annotation Reference
#    "Arc-Me"                         "Bladder"                       
#    "Bone-Marrow"                    "Bone-Marrow_c-kit"             
#    "Bone_Marrow_Mesenchyme"         "Brain"                         
#    "Bergman gliaBrain"              "E18-Brain"                     
#    "Embryonic-Mesenchyme"           "Embryonic-Stem-Cell"           
#    "Female_Fetal_Gonad"             "Fetal_Brain"                   
#    "Fetal_Intestine"                "Fetal_Kidney"                  
#    "Fetal-liver"                    "Fetal_Lung"                    
#    "Fetal_Stomache"                 "Kidney"                        
#    "Liver"                          "Lung"                          
#    "Lung-Mesenchyme"                "Male_Fetal_Gonad"              
#    "Mammary-Gland-Involution"       "Mammary-Gland-Lactation"       
#    "Mammary-Gland-Pregrancy"        "Mammary-Gland-Virgin"          
#    "Mesenchymal-Stem-Cell-Cultured" "Muscle"                        
#    "Neonatal_Brain"                 "Neonatal-Calvaria"             
#    "Neonatal-Heart"                 "Neonatal-Muscle"               
#    "Neonatal-Rib"                   "Neonatal-Skin"                 
#    "Ovary"                          "Pancreas"                      
#    "Peripheral_Blood"               "Placenta"                      
#    "Preimplantation-Embryo"         "Prostate"                      
#    "Retina"                         "Small-Intestine"               
#    "Spleen"                         "Stomach"                       
#    "Testis"                         "Thymus"                        
#    "Trophoblast-Stem-Cell"          "Uterus" 
