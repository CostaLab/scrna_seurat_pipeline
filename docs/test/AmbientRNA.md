---
sort: 3
---

# Ambient RNA

The pipeline has two possibilities to detect the amount of ambient RNA, SoupX and DecontX.
Right now, only DecontX is implemented in the master branch. For the SoupX, we need a different test data set. 
SoupX needs the filtered and unfiltered CellRanger results. DecontX only needs the Seruat object. 

## DecontX
DecontX needs the batch information per sample and determines the clustering per samples independently from the pipeline. 
The contamination and DEcontX clustering are then added as meta data to the Seurat object. 
A new assay with the decontamintated counts is added to the Seurat object, but is not used for the rest of the pipeline. 
We do not use the decontaminated reads for the pipeline, since the decontamination is not necessary reliable. 
