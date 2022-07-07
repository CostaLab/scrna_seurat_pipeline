library(tidyverse)

#https://panglaodb.se/markers.html?cell_type=%27all_cells%27
old_external <- read.csv("../external/Human_and_mouse_cell_markers-Markers.tsv", sep="\t")
all_mappings <- read.csv(file="../external/Human2Mouse_mappings.tsv", sep="\t")

dff <- read.csv("PanglaoDB_markers_27_Mar_2020.tsv", sep="\t")
dff <- dff %>% rename(Human.Gene=official.gene.symbol) %>%
        rename(Cell.Type =cell.type) %>%
        rename(Tissue.of.Origin=organ) %>%
        rename(Gene.Type=gene.type) %>%
        mutate(Added.By="PanglaoDB") %>%
        mutate(Added.Date="2022-07-05") %>%
        mutate(Reference="PanglaoDB") %>%
        mutate(Mouse.Gene="") %>%
        mutate(Comments="") %>%
        dplyr::filter(species=="Mm Hs")


dff <- dff[, colnames(old_external)]



### 1->multi, keep the title format if exists
all_to_del <- c()
for(Hs in unique(all_mappings$HGNC.symbol)){
  idxs <- which(all_mappings$HGNC.symbol == Hs)
  if(length(idxs) == 1){
    next
  }
  true_idx <- 1
  if(stringr::str_to_title(Hs) %in% all_mappings[idxs, "MGI.symbol"]){
    true_idx <- which(all_mappings[idxs, "MGI.symbol"] %in%  stringr::str_to_title(Hs))
  }
  del_idx <- idxs[-true_idx]

  all_to_del <- c(all_to_del,  del_idx)

}
keep <- setdiff(1:nrow(all_mappings), all_to_del)

u_all_mappings <- all_mappings[keep, ]


dff <- dff %>% dplyr::filter(Human.Gene %in% u_all_mappings$HGNC.symbol)
u_all_mappings <- u_all_mappings[u_all_mappings$HGNC.symbol %in% dff$Human.Gene, ]

dff$Mouse.Gene <- plyr::mapvalues(dff$Human.Gene,
                                  from=u_all_mappings$HGNC.symbol,
                                  to=u_all_mappings$MGI.symbol)

## delete same record if exists in old external
keys = paste0(dff$Human.Gene, "_", dff$Tissue.of.Origin)
dels <- which(keys %in% paste0(old_external$Human.Gene, "_", old_external$Tissue.of.Origin))
dff <- dff[-dels, ]


## remove useless
dff <- dff %>% dplyr::filter(!is.na(Tissue.of.Origin))
dff <- dff %>% dplyr::filter(length(Tissue.of.Origin)>0)


final_dff <- rbind(old_external, dff)

write.table(final_dff,  sep="\t", quote=F, file="panglaodb.tsv", row.names=F)



