library(scHCL)
library(Seurat)
library(stringr)
source("R/save_load_helper.R")



scrna <- load_object("exp/save/scrna_phase_clustering.Rds")
hcl_result <- scHCL(GetAssayData(object=scrna, slot="counts"), numbers_plot = 3)
corr=hcl_result$cors_matrix
rnms <- gsub("_[a-zA-Z]+\\.$", "", rownames(corr))
rnms <- gsub("_.*high", "", rnms)
rnms <- gsub("[0-9]\\.$", "", rnms)
rnms <- gsub("\\.[0-9]$", "", rnms)
rnms <- gsub("[0-9]\\.$", "", rnms)
rnms <- gsub("\\.\\.", ".", rnms)
rnms <- gsub("\\.$", "", rnms)

adults <- na.omit(unique(str_extract(rnms, "Adult.*$")))
fetal <- na.omit(unique(str_extract(rnms, "Fetal.*$")))

##adults
# [1] "Adult.Adipose"            "Adult.Adrenal.Gland"
# [3] "Adult.Artery"             "Adult.Ascending.Colon"
# [5] "Adult.Bladder"            "Adult.Bone.Marrow"
# [7] "Adult.Bone.Marrow.CD34N." "Adult.Bone.Marrow.CD34P."
# [9] "Adult.Brain"              "Adult.Cerebellum"
#[11] "Adult.Cervix"             "Adult.Duodenum"
#[13] "Adult.Epityphlon."        "Adult.Esophagus"
#[15] "Adult.Fallopian.Tube"     "Adult.Gall.Bladder"
#[17] "Adult.Heart"              "Adult.Ileum"
#[19] "Adult.JeJunum"            "Adult.Kidney"
#[21] "Adult.Liver"              "Adult.Lung"
#[23] "Adult.Muscle"             "Adult.Omentum"
#[25] "Adult.Pancreas"           "Adult.Peripheral.Blood"
#[27] "Adult.Pleura"             "Adult.Prostate"
#[29] "Adult.Rectum"             "Adult.Sigmoid.Colon"
#[31] "Adult.Spleen"             "Adult.Stomach"
#[33] "Adult.Temporal.Lobe"      "Adult.Thyroid"
#[35] "Adult.Trachea"            "Adult.Transverse.Colon"
#[37] "Adult.Ureter"             "Adult.Uterus"
#[39] "Adult.hepatocyte.Liver"
##Fetal
# [1] "Fetal.Adrenal.Gland"                "Fetal.Brain"
# [3] "Fetal.Calvaria"                     "Fetal.Eyes"
# [5] "Fetal.Female.Gonad"                 "Fetal.Gonad"
# [7] "Fetal.Heart"                        "Fetal.Intestine"
# [9] "Fetal.Kidney"                       "Fetal.Liver"
#[11] "Fetal.Lung"                         "Fetal.Male.Gonad"
#[13] "Fetal.leydig.cell.Fetal.Male.Gonad" "Fetal.Mid.Brain"
#[15] "Fetal.Muscle"                       "Fetal.Pancreas"
#[17] "Fetal.Rib"                          "Fetal.Skin"
#[19] "Fetal.Spinal.Cord"                  "Fetal.Stomach"
#[21] "Fetal.Thymus"                       "Fetal.hepatocyte.Liver"



idx1 <- which(grepl("Adult", rnms))
idx2 <- which(grepl("Fetal", rnms))
rrnms <- rnms[c(-idx1, -idx2)]
rrnms <- na.omit(unique(gsub("^\\.cell\\.", "", str_extract(rrnms, "\\.cell.*$"))))
rrnms <- unique(gsub("^\\.cell[s]*[0-9]*\\.", "", rrnms))

## Others
# [1] "Airway.Epithelium"              "Breast.Epithelium"
# [3] "Chorionic.Villus"               "Cord.Blood"
# [5] "Centrocyte.Cord.Blood"          "Unknown.Cord.Blood"
# [7] "Plasmocyte.Cord.Blood"          "Cord.Blood.CD34P"
# [9] "Dendritic.Cell.Monocyte"        "ES.to.EB_8Day"
#[11] "Haematopoietic.Stem.Cell"       "Liver"
#[13] "Neonatal.Adrenal.Gland"         "Placenta"
#[15] "embryo.Preimplantation.Embryo." "Testis"
