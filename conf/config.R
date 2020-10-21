### --------------Initail info----------------------------
PROJECT = "Mouse heart Gli&CD45 project" ## set project name 
ORGAN = 'Heart'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"         #For external annotation. Options: Human, Mouse
MCA_NAME = "Neonatal-Heart" #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/ 

# filtering params when create seurat object
MINCELLS  = 5 
MINGENES  = 50

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv" 

data_src = c( 

     NK1_Gli1_IRI     =   "data/summed_mtx/sum_NK1_Gli1_IRI/",
     NK2_CD45_IRI     =   "data/summed_mtx/sum_NK2_CD45_IRI/",
     NK3_Gli1_Sham    =   "data/summed_mtx/sum_NK3_Gli1_Sham/",
     NK4_CD45_Sham    =   "data/summed_mtx/sum_NK4_CD45_Sham/"
)


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
       scrna_mitoRegressOut       = 1, ## Regress out mito only
       scrna_riboRegressOut       = 0, ## Regress out ribo only
       scrna_mitoRiboRegressOut   = 0, ## Regress out ribo&mito 
       scrna_RegressOutAll        = 0, ## Regress out celly cycle & mito & ribo
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
