
# intro -------------------------------------------------------------------

# immune profiling with scRNAseq data
# revision to NP137 manuscript

# cell:cell communication on scRNAseq data



# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")

suppressPackageStartupMessages({
  library(Seurat);  library(SeuratDisk);#; library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr); library(xlsx);  library(stringr)
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(biomaRt)
  library(Scillus)
  library(SCENIC)
  library(snow)
  library(DoubletFinder)
  library(Nebulosa)
  library(celldex)
  library(ProjecTILs)
  library(ComplexHeatmap)
  
  
})

set.seed(170323)


# CellChat grouped -----------------------------------

rm(list=ls())

load(file="sc.all.celltypes.RData")

DimPlot(sc.merged)

meta.Tcells <- read.csv(file="metadata.Tcells.csv", row.names = 1)
head(meta.Tcells)

#load("epithelial.RData")
#load("Tcells.RData")
#DimPlot(epithelial)
#DimPlot(Tcells)

Idents(epithelial) <- "epithelial"
epithelial$cellchat.cat <- Idents(epithelial)
Idents(Tcells)
Tcells$cellchat.cat <- Idents(Tcells)
DefaultAssay(Tcells) <- "SCT"

sc.merged <- merge(epithelial, Tcells)
levels(sc.merged)
head(sc.merged[[]])

#

Tex <- rownames(meta.Tcells[meta.Tcells$final.cats == "CD8 Tex", ])
Tregs <- rownames(meta.Tcells[meta.Tcells$final.cats == "Treg cells", ])
CD8 <- rownames(meta.Tcells[meta.Tcells$final.cats == "CD8 T cells", ])
CD4 <- rownames(meta.Tcells[meta.Tcells$final.cats == "CD4 T cells", ])
NK <- rownames(meta.Tcells[meta.Tcells$final.cats == "Natural killer cells", ])

sc.merged <- SetIdent(sc.merged, cells = Tex, value = 'CD8 Tex')
sc.merged <- SetIdent(sc.merged, cells = Tregs, value = 'Treg cells')
sc.merged <- SetIdent(sc.merged, cells = CD8, value = 'CD8 T cells')
sc.merged <- SetIdent(sc.merged, cells = CD4, value = 'CD4 T cells')
sc.merged <- SetIdent(sc.merged, cells = NK, value = 'NK cells')

DimPlot(sc.merged)
sc.merged <- subset(sc.merged, idents = "CD8 NKT-like cells", invert = T)
sc.merged$cellchat.cat <- Idents(sc.merged)
table(sc.merged$cellchat.cat)

library(CellChat)
cellchat <- createCellChat(object = sc.merged, 
                           group.by = "cellchat.cat")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 100)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

jpeg("netVisual_circle_all.cells.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  jpeg(paste0("netVisual_circle", rownames(mat)[i], ".jpeg"), height = 900, width = 900, quality = 100)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

jpeg("netAnalysis_heatmap_all.cells.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat, file = "cellchat_all.cells.RData")
save(sc.merged, file="sc.merged.cellchat.RData")


# individual 1 ----------------------------------

rm(list=ls())

load("sc.merged.cellchat.RData")

NP3 <- subset(sc.merged, sample == "NP3")
NP4 <- subset(sc.merged, sample == "NP4")


library(CellChat)
cellchat.1 <- createCellChat(object = NP3, group.by = "cellchat.cat")
cellchat.1@DB <- CellChatDB.human
cellchat.1 <- subsetData(cellchat.1)
future::plan("multiprocess", workers = 4)
cellchat.1 <- identifyOverExpressedGenes(cellchat.1)
cellchat.1 <- identifyOverExpressedInteractions(cellchat.1)
cellchat.1 <- projectData(cellchat.1, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.1 <- computeCommunProb(cellchat.1)
cellchat.1 <- filterCommunication(cellchat.1, min.cells = 100)
cellchat.1 <- computeCommunProbPathway(cellchat.1)
cellchat.1 <- aggregateNet(cellchat.1)
groupSize <- as.numeric(table(cellchat.1@idents))
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", 
                                  font.size = 9, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", 
                                    font.size = 9, width = 4, height = 10)

jpeg("NP3_netVisual_circle_all.cells.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
jpeg("NP3_netAnalysis_heatmap_all.cells.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.1, file = "cellchat.NP3.RData")


cellchat.2 <- createCellChat(object = NP4, group.by = "cellchat.cat")
cellchat.2@DB <- CellChatDB.human
cellchat.2 <- subsetData(cellchat.2)
future::plan("multiprocess", workers = 4)
cellchat.2 <- identifyOverExpressedGenes(cellchat.2)
cellchat.2 <- identifyOverExpressedInteractions(cellchat.2)
cellchat.2 <- projectData(cellchat.2, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.2 <- computeCommunProb(cellchat.2)
cellchat.2 <- filterCommunication(cellchat.2, min.cells = 100)
cellchat.2 <- computeCommunProbPathway(cellchat.2)
cellchat.2 <- aggregateNet(cellchat.2)
groupSize <- as.numeric(table(cellchat.2@idents))
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", 
                                  font.size = 8, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", 
                                    font.size = 8, width = 4, height = 10)

jpeg("NP4_netVisual_circle_all.cells.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
jpeg("NP4_netAnalysis_heatmap_all.cells.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.2, file = "cellchat.NP4.RData")





# individual 2 ----------------------------------

rm(list=ls())

load("sc.merged.cellchat.RData")

table(Idents(sc.merged))

# select 
sc.merged <- subset(sc.merged, idents = c("NK cells", "CD4 T cells", "CD8 T cells",
                                          "Treg cells", "Tumor cells", "Neutrophils",
                                          "B cells", "Macrophages", "Endothelial",
                                          "CAFs"), invert = F)
sc.merged$cellchat.cat <- Idents(sc.merged)

NP3 <- subset(sc.merged, sample == "NP3")
NP4 <- subset(sc.merged, sample == "NP4")


library(CellChat)
cellchat.1 <- createCellChat(object = NP3, group.by = "cellchat.cat")
cellchat.1@DB <- CellChatDB.human
cellchat.1 <- subsetData(cellchat.1)
future::plan("multiprocess", workers = 4)
cellchat.1 <- identifyOverExpressedGenes(cellchat.1)
cellchat.1 <- identifyOverExpressedInteractions(cellchat.1)
cellchat.1 <- projectData(cellchat.1, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.1 <- computeCommunProb(cellchat.1)
cellchat.1 <- filterCommunication(cellchat.1, min.cells = 100)
cellchat.1 <- computeCommunProbPathway(cellchat.1)
cellchat.1 <- aggregateNet(cellchat.1)
groupSize <- as.numeric(table(cellchat.1@idents))
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", 
                                  font.size = 9, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", 
                                    font.size = 9, width = 4, height = 10)

jpeg("NP3_netVisual_circle_2.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
jpeg("NP3_netAnalysis_heatmap_2.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.1, file = "cellchat.NP3_2.RData")


cellchat.2 <- createCellChat(object = NP4, group.by = "cellchat.cat")
cellchat.2@DB <- CellChatDB.human
cellchat.2 <- subsetData(cellchat.2)
future::plan("multiprocess", workers = 4)
cellchat.2 <- identifyOverExpressedGenes(cellchat.2)
cellchat.2 <- identifyOverExpressedInteractions(cellchat.2)
cellchat.2 <- projectData(cellchat.2, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.2 <- computeCommunProb(cellchat.2)
cellchat.2 <- filterCommunication(cellchat.2, min.cells = 100)
cellchat.2 <- computeCommunProbPathway(cellchat.2)
cellchat.2 <- aggregateNet(cellchat.2)
groupSize <- as.numeric(table(cellchat.2@idents))
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", 
                                  font.size = 8, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", 
                                    font.size = 8, width = 4, height = 10)

jpeg("NP4_netVisual_circle_2.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
jpeg("NP4_netAnalysis_heatmap_2.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.2, file = "cellchat.NP4_2.RData")





# individual 3 ----------------------------------

rm(list=ls())

load("sc.merged.cellchat.RData")

table(Idents(sc.merged))

# select 
sc.merged <- subset(sc.merged, idents = c("NK cells", "CD4 T cells", "CD8 T cells",
                                          "Treg cells", "Tumor cells"), invert = F)
sc.merged$cellchat.cat <- Idents(sc.merged)

NP3 <- subset(sc.merged, sample == "NP3")
NP4 <- subset(sc.merged, sample == "NP4")


library(CellChat)
cellchat.1 <- createCellChat(object = NP3, group.by = "cellchat.cat")
cellchat.1@DB <- CellChatDB.human
cellchat.1 <- subsetData(cellchat.1)
future::plan("multiprocess", workers = 8)
cellchat.1 <- identifyOverExpressedGenes(cellchat.1)
cellchat.1 <- identifyOverExpressedInteractions(cellchat.1)
cellchat.1 <- projectData(cellchat.1, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.1 <- computeCommunProb(cellchat.1)
cellchat.1 <- filterCommunication(cellchat.1, min.cells = 100)
cellchat.1 <- computeCommunProbPathway(cellchat.1)
cellchat.1 <- aggregateNet(cellchat.1)
groupSize <- as.numeric(table(cellchat.1@idents))
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", 
                                  font.size = 9, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", 
                                    font.size = 9, width = 4, height = 10)

jpeg("NP3_netVisual_circle_3.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")
jpeg("NP3_netAnalysis_heatmap_3.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.1, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.1, file = "cellchat.NP3_3.RData")


cellchat.2 <- createCellChat(object = NP4, group.by = "cellchat.cat")
cellchat.2@DB <- CellChatDB.human
cellchat.2 <- subsetData(cellchat.2)
future::plan("multiprocess", workers = 8)
cellchat.2 <- identifyOverExpressedGenes(cellchat.2)
cellchat.2 <- identifyOverExpressedInteractions(cellchat.2)
cellchat.2 <- projectData(cellchat.2, PPI.human)
options(future.globals.maxSize = 8000 * 1024^2)
cellchat.2 <- computeCommunProb(cellchat.2)
cellchat.2 <- filterCommunication(cellchat.2, min.cells = 100)
cellchat.2 <- computeCommunProbPathway(cellchat.2)
cellchat.2 <- aggregateNet(cellchat.2)
groupSize <- as.numeric(table(cellchat.2@idents))
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", 
                                  font.size = 8, width = 4, height = 10) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", 
                                    font.size = 8, width = 4, height = 10)

jpeg("NP4_netVisual_circle_3.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
jpeg("NP4_netAnalysis_heatmap_3.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.2, file = "cellchat.NP4_3.RData")




# get tables ---------------------------------------------------------

rm(list=ls())

# function to load .RData objects on a different name
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

base.dir <- "cellchat/individual_3/"
files <- list.files(path= base.dir, pattern = "RData")
names <- str_sub(files, 10, -9)

for(i in 1:length(files)){
  cellchat <- loadRData(paste0(base.dir, files[i]))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  df.net <- subsetCommunication(cellchat)
  df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
  write.csv(df.net, file = paste0(base.dir, "net_", names[i] , ".csv"))
  write.csv(df.net.2, file = paste0(base.dir, "net2_", names[i] , ".csv"))
  saveRDS(cellchat, file = paste0(base.dir, "cellchat_", names[i], ".rds"))
}


# compare tables

sc1 <- read.csv("cellchat/individual_3/net2_NP3.csv", row.names = 1)
table(sc1$source)
sc1 <- sc1[sc1$source == "Tumor cells", ]
sc1 <- sc1[sc1$target == "CD8 T cells", ]
#sc1 <- sc1[sc1$target == "Tumor cells", ]
#sc1 <- sc1[sc1$source == "CD8 T cells", ]
head(sc1)

sc2 <- read.csv("cellchat/individual_3/net2_NP4.csv", row.names = 1)
table(sc2$source)
sc2 <- sc2[sc2$source == "Tumor cells", ]
sc2 <- sc2[sc2$target == "CD8 T cells", ]
#sc2 <- sc2[sc2$target == "Tumor cells", ]
#sc2 <- sc2[sc2$source == "CD8 T cells", ]
head(sc2)

intersect(sc1$pathway_name, sc2$pathway_name)
setdiff(sc1$pathway_name, sc2$pathway_name)
setdiff(sc2$pathway_name, sc1$pathway_name)


sp1a <- read.csv("spatial/cellchat.2/net2_P034_Pre.csv", row.names = 1)
table(sp1a$source)
sp1a <- sp1a[sp1a$source == "Tumor cells", ]
sp1a <- sp1a[sp1a$target == "T cells", ]
head(sp1a)

sp2a <- read.csv("spatial/cellchat.2/net2_P034_Post.csv", row.names = 1)
sp2a <- sp2a[sp2a$source == "Tumor cells", ]
sp2a <- sp2a[sp2a$target == "T cells", ]

intersect(sc1$pathway_name, sp1a$pathway_name)

sp1b <- read.csv("spatial/cellchat.2/net2_P039_Pre.csv", row.names = 1)
table(sp1b$source)
sp1b <- sp1b[sp1b$source == "Tumor cells", ]
sp1b <- sp1b[sp1b$target == "T cells", ]
head(sp1b)

sp2b <- read.csv("spatial/cellchat.2/net2_P039_Post.csv", row.names = 1)
sp2b <- sp2b[sp2b$source == "Tumor cells", ]
sp2b <- sp2b[sp2b$target == "T cells", ]

intersect(sc2$pathway_name, sp1b$pathway_name)
intersect(sp1a$pathway_name, sp1b$pathway_name)
setdiff(sp1a$pathway_name, sp1b$pathway_name)
setdiff(sp1b$pathway_name, sp1a$pathway_name)

int1 <- intersect(sp2b$pathway_name, 
          intersect(sc2$pathway_name, sp2a$pathway_name))
# "COLLAGEN" "LAMININ"  "MHC-I"    "MIF"      "MK" 

sc1[sc1$pathway_name%in%int1, ]
sc2[sc2$pathway_name%in%int1, ]
sp1a[sp1a$pathway_name%in%int1, ]
sp2a[sp2a$pathway_name%in%int1, ]
sp1b[sp1b$pathway_name%in%int1, ]
sp2b[sp2b$pathway_name%in%int1, ]

# MHC-I seems the most consistent finding


# comparative ----------------------------------

rm(list=ls())

cellchat.1 <- readRDS("cellchat/individual_3/cellchat_NP3.rds")
cellchat.2 <- readRDS("cellchat/individual_3/cellchat_NP4.rds")

object.list <- list('NP3' = cellchat.1, 'NP4' = cellchat.2)

# Compare the number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count",
                           size.text = 14,
                           title.name = "Number of Interactions")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",
                           size.text = 14,
                           title.name = "Interaction Strength")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(1, 2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 2))
# The width of edges represent the relative number of interactions or interaction strength. 
# Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(1, 2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1, 2))
gg1 + gg2

# To better control the node size and edge weights of the inferred networks across different datasets, 
# we compute the maximum number of cells per cell group and the maximum number of interactions 
# (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


# Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
levels(cellchat@idents$joint)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor cells", comparison = c(1, 2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8 T cells", comparison = c(1, 2))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 T cells", comparison = c(1, 2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK cells", comparison = c(1, 2))
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Treg cells", comparison = c(1, 2))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4,gg5))


# Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", min_dist = 0.5) # Sensible values are in the range 0.001 to 0.5
cellchat <- netClustering(cellchat, type = "functional", k=10)
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, 
                            type = "functional", 
                            dot.size = c(10, 10), label.size = 5, 
                            title = "Functional Similarity")

rankSimilarity(cellchat, type = "functional")

# Compare the overall information flow of each signaling pathway
# 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2), font.size = 8)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2), font.size = 8)
gg1 + gg2


# Compare outgoing (or incoming) signaling

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 10)
draw(ht1 + ht2, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(2, "cm"))


# Identify the upregulated and down-regulated signaling ligand-receptor pairs

netVisual_bubble(cellchat, 
                 #                 sources.use = 1, 
                 #                 targets.use = 1:3,  
                 comparison = c(1, 2), 
                 font.size = 10,
                 thresh = 0.05, remove.isolate = T,
                 angle.x = 45)

levels(cellchat@idents$joint)

gg1 <- netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 5,  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in NP4", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 5,  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in NP4", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4),  
                        comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in NP4", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4),  
                        comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in NP4", angle.x = 45, remove.isolate = T)
gg1 + gg2


# Identify dysfunctional signaling by using differential expression analysis

# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups.
# Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. 
# Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, 
# and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells. Such analysis can be done as follows.

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "NP4"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "NP4",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NP3",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5, targets.use = c(1:4), font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5, targets.use = c(1:4), font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:4), targets.use = 5, font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:4), targets.use = 5, font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = c(1:4), 
                     targets.use = 5, 
                     slot.name = 'net', 
                     net = net.up, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], 
                     sources.use = 1:4, 
                     targets.use = 5, 
                     slot.name = 'net', 
                     net = net.down, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = 5, 
                     targets.use = 1:4, 
                     slot.name = 'net', 
                     net = net.up, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], 
                     sources.use = 5, 
                     targets.use = 1:4, 
                     slot.name = 'net', 
                     net = net.down, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human', measure = "signaling")

# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', measure = "signaling")

# Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NP3", "NP4")) # set factor level
plotGeneExpression(cellchat, signaling = "MIF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "NECTIN", split.by = "datasets", colors.ggplot = T)



# global communication patterns

# river plot
library(ggalluvial)
# selectK(cellchat.1, pattern = "outgoing") # takes too long
par(mfrow=c(1,1))

cellchat.1 <- identifyCommunicationPatterns(cellchat.1, pattern = "outgoing", k = 4)
netAnalysis_river(cellchat.1, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

cellchat.2 <- identifyCommunicationPatterns(cellchat.2, 
                                            pattern = "outgoing", 
                                            slot.name = "netP",
                                            k = 3)
netAnalysis_river(cellchat.2, pattern = "outgoing")
netAnalysis_dot(cellchat.2, pattern = "outgoing")

cellchat.2 <- identifyCommunicationPatterns(cellchat.2, 
                                            pattern = "incoming", 
                                            slot.name = "netP",
                                            k = 3)
netAnalysis_river(cellchat.2, pattern = "incoming")
netAnalysis_dot(cellchat.2, pattern = "incoming")


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}





df.net <- subsetCommunication(cellchat)
head(df.net)
df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
head(df.net.2)

write.csv(df.net$NP3, file = "cellchat.net.1.NP3.csv")
write.csv(df.net$NP4, file = "cellchat.net.1.NP4.csv")
write.csv(df.net.2$NP3, file = "cellchat.net.2.NP3.csv")
write.csv(df.net.2$NP4, file = "cellchat.net.2.NP4.csv")


# other viz

#pathways.show <- c("TGFb") 
pathways.show <- c("COLLAGEN") 

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to TAMs and the right portion shows signaling to all other cells 
levels(cellchat@idents)
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver, 
                    layout = "hierarchy")

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1:17, targets.use = c(14:15), remove.isolate = F)
netVisual_bubble(cellchat, sources.use = 14:15, targets.use = 1:17, remove.isolate = F)

# network centrality score for a given pathway
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 15, height = 6, font.size = 12)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MIF"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together


selectK(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")

nPatterns = 4
#nPatterns = 8
par(mar=c(10,5,5,5), mfrow=c(1,1))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", 
                                          width = 10, height = 20, font.size = 10,
                                          k = nPatterns)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", 
                                          width = 10, height = 20, font.size = 10,
                                          k = nPatterns)

# river plot
netAnalysis_river(cellchat.1, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "incoming")






# end ---------------------------------------------------------------------







