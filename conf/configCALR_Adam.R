### --------------Initial info----------------------------
PROJECT = "CALR_Adam" ## set project name
ORGAN = 'Blood'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Bone-Marrow" #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/
HCL_NAME = "Adult-Bone-Marrow-CD34P" #For HCL annotation.


# filtering params when create seurat object
MINCELLS  = 5
MINGENES  = 50


INTEGRATION_OPTION = "seurat" ### or harmony

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv"

data_src = c(
      "Empty EV #31-1"      = "data/A_MxCre",
      "Empty EV #31-2"      = "data/B_MxCre",
      "CALR deletion #34-1" = "data/C_Csnk",
      "CALR deletion #34-2" = "data/D_Csnk"
)



##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        "Empty EV #31-1" = "Empty EV",
	"Empty EV #31-2" = "Empty EV",
        "CALR deletion #34-1" = "CALR deletion",
        "CALR deletion #34-2" = "CALR deletion"
)



## Phase_1, set 1 to regressout
preprocess_regressout = c("mito"       = 1,
                          "ribo"       = 0,
                          "cellcycle"  = 1,
			  "ambientRNA" = 1)


#Analysis_phases
#1. scrna_phase_preprocess
#2. scrna_phase_clustering
#3. scrna_phase_comparing


### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load
conf = c(
         scrna_phase_clustering     = 2, ## integration & clustering
	 scrna_ambient_rna          = 1, ## estimating and correcting for ambient RNA
         scrna_remove_clusters      = 1, ## remove clusters
         scrna_cluster_annotation   = 1, ## Annotate clusters according to `cluster_annotation`
         scrna_phase_comparing      = 1 ## DE GO pathway analysis etc. All rest calculating will be stored here
        )

### ----------specific settings for some functions ----------------

scrna_remove_clusters = c(3) # Cluster 3 is the erythrocytes cluster


### cluster annotation
from_cluster_slot = "removed_clusters"
cluster_annotation <- c(
    "0" = "MSC-2",
    "1" = "MSC",
    "2" = "MSC-1",
    "4" = "MSC-3",
    "5" = "SCP",
    "6" = "OLC",
    "7" = "Unknown1",
    "8" = "Unknown2"
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
