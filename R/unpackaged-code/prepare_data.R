
######################################

human_TFs <- read.csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt", header = F)
human_TFs <- human_TFs$V1
print(length(human_TFs))
human.tfs <- as.character(human_TFs)

human_TF2 <- read.csv("/Users/zxli/Downloads/_TF.txt", sep = "\t")
human.tfs <- unique(c(human.tfs, as.character(human_TF2$Symbol)))
print(length(human.tfs))

save(human.tfs, file = "/Users/zxli/Dropbox/bigWriting/8.iterbi/R_package/iterbi/data/human.tfs.rda")

######################################

mouse.tfs <- read.csv("/Users/zxli/Downloads/Mus_musculus_TF.txt", sep = "\t")$Symbol
mouse.tfs <- as.character(mouse.tfs)
print(length(mouse.tfs))

save(mouse.tfs, file = "/Users/zxli/Dropbox/bigWriting/8.iterbi/R_package/iterbi/data/mouse.tfs.rda")


######################################
library(CellChat)
head(CellChatDB.human$interaction)

tmp <- CellChatDB.human$interaction$interaction_name_2

ligand <- unlist(lapply(tmp, function(x){
  strsplit(x, split = " - ")[[1]][1]
} ))
human.ligand <- unique(ligand)
human.ligand <- unique(unlist(strsplit(human.ligand, split = "_")))

receptors <- unlist(lapply(tmp, function(x){
  strsplit(x, split = " - ")[[1]][2]
} ))
receptors_1 <- unique(receptors[!grepl("\\+", receptors)])
receptors_2 <- unique(unlist(strsplit(receptors[grepl("\\+", receptors)], split = "\\(|\\+|\\)")))
receptors_2 <- receptors_2[!receptors_2 %in% ""]
human.receptor <- unique(c(receptors_1, receptors_2))
length(human.receptor)

human.ligand <- trimws(human.ligand)
human.receptor <- trimws(human.receptor)
human.ligand.receptor <- list(`human.ligand`=human.ligand, `human.receptor`=human.receptor)

save(human.ligand.receptor, file = "/Users/zxli/Dropbox/bigWriting/8.iterbi/R_package/iterbi/data/human.ligand.receptor.rda")

######################################
tmp <- CellChatDB.mouse$interaction$interaction_name_2

ligand <- unlist(lapply(tmp, function(x){
  strsplit(x, split = " - ")[[1]][1]
} ))
mouse.ligand <- unique(ligand)
mouse.ligand <- unique(unlist(strsplit(mouse.ligand, split = "_")))

receptors <- unlist(lapply(tmp, function(x){
  strsplit(x, split = " - ")[[1]][2]
} ))
receptors_1 <- unique(receptors[!grepl("\\+", receptors)])
receptors_2 <- unique(unlist(strsplit(receptors[grepl("\\+", receptors)], split = "\\(|\\+|\\)")))
receptors_2 <- receptors_2[!receptors_2 %in% ""]
mouse.receptor <- unique(c(receptors_1, receptors_2))
length(mouse.receptor)

mouse.ligand <- trimws(mouse.ligand)
mouse.receptor <- trimws(mouse.receptor)
mouse.ligand.receptor <- list(`mouse.ligand`=mouse.ligand, `mouse.receptor`=mouse.receptor)

save(mouse.ligand.receptor, file = "/Users/zxli/Dropbox/bigWriting/8.iterbi/R_package/iterbi/data/mouse.ligand.receptor.rda")
