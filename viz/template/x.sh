#sed -i "s/scrna <- readRDS.*$//" DE-GO-1v1.template
#sed -i "s/scrna <- readRDS.*$//"  DE-GO-stagesVS.template
#sed -i "s/scrna <- readRDS.*$//"  Genesets-1v1.template
#sed -i "s/scrna <- readRDS.*$//"  Genesets-stageVS.template
#sed -i "s/scrna <- readRDS.*$//"  hallmark-1v1.template
#sed -i "s/scrna <- readRDS.*$//"  hallmark-stageVS.template
#sed -i "s/scrna <- readRDS.*$//"  kegg-1v1.template
#sed -i "s/scrna <- readRDS.*$//"  kegg-stageVS.template
#sed -i "s/scrna <- readRDS.*$//"  progeny-stageVS.template
#sed -i "s/scrna <- readRDS.*$//"  reactome-1v1.template
#sed -i "s/scrna <- readRDS.*$//"  reactome-stageVS.template


tmpls=(DE-GO-1v1.template
       DE-GO-stagesVS.template
       Genesets-1v1.template
       Genesets-stageVS.template
       hallmark-1v1.template
       hallmark-stageVS.template
       index.template
       kegg-1v1.template
       kegg-stageVS.template
       progeny-stageVS.template
       reactome-1v1.template
       reactome-stageVS.template
)


for f in "${tmpls[@]}"; do
     sed -i 's/library(ggplot2)/suppressPackageStartupMessages(library(ggplot2))/' $f
     sed -i 's/library(gridExtra)/suppressPackageStartupMessages(library(gridExtra))/' $f
     sed -i 's/library(cowplot)/suppressPackageStartupMessages(library(cowplot))/' $f
     sed -i 's/library(dplyr)/suppressPackageStartupMessages(library(dplyr))/' $f
     sed -i 's/library(openxlsx)/suppressPackageStartupMessages(library(openxlsx))/' $f
     sed -i 's/library(WriteXLS)/suppressPackageStartupMessages(library(WriteXLS))/' $f
     sed -i 's/library(stringr)/suppressPackageStartupMessages(library(stringr))/' $f
     sed -i 's/library(digest)/suppressPackageStartupMessages(library(digest))/' $f
     sed -i 's/library(ggridges)/suppressPackageStartupMessages(library(ggridges))/' $f
     sed -i 's/library(plotly)/suppressPackageStartupMessages(library(plotly))/' $f
     sed -i 's/library(Seurat)/suppressPackageStartupMessages(library(Seurat))/' $f
     sed -i 's/library(Hmisc)/suppressPackageStartupMessages(library(Hmisc))/' $f
     sed -i 's/library(EnhancedVolcano)/suppressPackageStartupMessages(library(EnhancedVolcano))/' $f
     sed -i 's/library(ComplexHeatmap)/suppressPackageStartupMessages(library(ComplexHeatmap))/' $f
     sed -i 's/library(glue)/suppressPackageStartupMessages(library(glue))/' $f
done
