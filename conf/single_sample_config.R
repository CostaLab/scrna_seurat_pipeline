### --------------Initail info----------------------------
PROJECT = "Mouse pbmc project" ## set project name 
ORGAN = 'Blood'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Bone-Marrow"  #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/ 

# filtering params when create seurat object
MINCELLS  = 5 
MINGENES  = 50

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv" 

data_src = c( 
     pbmc     =   "data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
)

##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
    pbmc      = "pbmc"
)

### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load   3.  
conf = c(
       scrna_rawdata              = 1, ## read count matrix and merge samples to a Seurat obj 
       scrna_filter               = 1, ## filter nFeatureRNA and nCountRNA
       scrna_preprocess           = 1, ## Normalize & FindvariablegFeatures and ScaleData
       scrna_cellcycle            = 1, ## Cell cycle scoring
       scrna_cycleRegressOut      = 1, ## Regress out cell cycle effects
       scrna_RegressOutAll        = 1, ## Regress out celly cycle & mito & ribo
       scrna_sltn_batch_clustering= 1, ## only one sample
       scrna_singleton_clustering = 1) ## clustering 




### ----------specific settings for some functions ----------------



scrna_clusterlevel_filter_df = do.call(rbind, list(
  mito_cluster1_filter  =  data.frame(cluster=1, max_pct=1, min_pct=0, type="mito"),
  mito_cluster2_filter  =  data.frame(cluster=1, max_pct=1, min_pct=0, type="ribo"),
  mito_cluster3_filter  =  data.frame(cluster=1, max_pct=1, min_pct=0, type="exclude"),
  mito_cluster4_filter  =  data.frame(cluster=1, max_pct=1, min_pct=0, type="mito"),
  mito_cluster5_filter  =  data.frame(cluster=1, max_pct=1, min_pct=0, type="mito")
)) 
                           
                                 




scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)


scrna_remove_clusters = c(1, 3, 6)
scrna_remove_recluster = c(1, 3, 6)



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
