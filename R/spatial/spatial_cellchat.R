
# spatial RNAseq data
# NP137 manuscript
# Visium for two samples, before and after therapy

# spatial cellchat


# libraries ---------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(Seurat)
  library(ComplexHeatmap)
  library(stringr)
  
})

options(stringsAsFactors = FALSE)

set.seed(4)

setwd()


# run cell chat in each sample --------------------------------------------------------------------

load("st.list.anno.manual.RData")
names <- names(new.anno.list)
image.dirs <- list.files(path = "/Visium", full.names = T)

# CD4 and CD8 T cells are grouped because of their low numbers

for(i in 1:length(new.anno.list)){
  
  sc.query <- new.anno.list[[i]]
  Idents(sc.query) <- sc.query$cell_type
  sc.query <- RenameIdents(sc.query, 'CD8 T cells' = 'T cells', 'CD4 T cells' = 'T cells')
  sc.query <- subset(sc.query, idents = 'Undefined', invert = T)
  
  # Prepare input data for CelChat analysis
  data.input = GetAssayData(sc.query, slot = "data", assay = "SCT") # normalized data matrix
  meta = data.frame(labels = Idents(sc.query), 
                    row.names = names(Idents(sc.query))) # manually create a dataframe consisting of the cell labels

  # load spatial imaging information
  spatial.locs = GetTissueCoordinates(sc.query, scale = NULL, cols = c("imagerow", "imagecol")) 
  scale.factors = jsonlite::fromJSON(txt = file.path(paste0(image.dirs[i], "/outs/spatial"), 'scalefactors_json.json'))
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )

  cellchat <- createCellChat(object = data.input, 
                             meta = meta, 
                             group.by = "labels",
                             datatype = "spatial", 
                             coordinates = spatial.locs, 
                             scale.factors = scale.factors)

  # set the used database in the object
  cellchat@DB <- CellChatDB.human
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 8) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI 
  # (Optional: when running it, USER should set `raw.use = FALSE` 
  # in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5) 
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))

  # visualization

  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, 
                   vertex.weight = rowSums(cellchat@net$count), 
                   weight.scale = T, label.edge= F, 
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, 
                   vertex.weight = rowSums(cellchat@net$weight), 
                   weight.scale = T, label.edge= F, 
                   title.name = "Interaction weights/strength")
  
  save(cellchat, file = paste0("cellchat_", names[i], ".RData"))

  rm(sc.query, data.input, meta, spatial.locs, scale.factors, groupSize, cellchat)
  
  
}




# compare by patient -------------------------------------------------------

rm(list=ls())

load("cellchat.2/cellchat_P034_Pre.RData")
#load("cellchat.2/cellchat_P039_Pre.RData")
cellchat.1 <- cellchat
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1, slot.name = "netP")

load("cellchat.2/cellchat_P034_Post.RData")
#load("cellchat.2/cellchat_P039_Post.RData")
cellchat.2 <- cellchat
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")

object.list <- list('Pre' = cellchat.1, 'Post' = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

# Compare the number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count",
                           size.text = 14,
                           title.name = "Number of Interactions")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",
                           size.text = 14,
                           title.name = "Interaction Strength")
gg1 + gg2
ggsave(filename = "comparative.number.of.interactions.jpeg")

#jpeg("differential.number.of.interactions.jpeg", quality = 100)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", comparison = c(1, 2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 2))
# The width of edges represent the relative number of interactions or interaction strength. 
# Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#dev.off()

gg1 <- netVisual_heatmap(cellchat, measure = "count", comparison = c(1, 2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1, 2))
gg1 + gg2
#ggsave(filename = "differential.heatmap.jpeg", plot = gg1)

# To better control the node size and edge weights of the inferred networks across different datasets, 
# we compute the maximum number of cells per cell group and the maximum number of interactions 
# (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
#jpeg("chord.number.of.interactions.jpeg", quality = 100)
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
ggsave("outgoing.incoming.PCA.jpeg")

levels(cellchat@idents$joint)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor cells", comparison = c(1, 2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T cells", comparison = c(1, 2))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial", comparison = c(1, 2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophages", comparison = c(1, 2))
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CAFs", comparison = c(1, 2))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4,gg5))
ggsave("outgoing.incoming.by.celltype.jpeg")


# Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", min_dist = 0.5) # Sensible values are in the range 0.001 to 0.5
cellchat <- netClustering(cellchat, type = "functional", k=10)
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, 
                            type = "functional", 
                            dot.size = c(10, 10), label.size = 5, 
                            title = "Functional Similarity")
ggsave("functional.similarity.jpeg")


rankSimilarity(cellchat, type = "functional", font.size = 12)
ggsave("ranked.similarity.jpeg")

# Compare the overall information flow of each signaling pathway
# 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2), font.size = 8)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2), font.size = 8)
gg1 + gg2
ggsave("compared.rankings.jpeg")
ggsave("compared.rankings.pdf")


# Compare outgoing (or incoming) signaling

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 25)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 25)
#jpeg("signals.heatmap.jpeg", quality = 100)
draw(ht1 + ht2, ht_gap = unit(2, "cm"))
#dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 25, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 25, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 25, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 25, color.heatmap = "OrRd")
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

gg1 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = 4,  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in NP4", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = 4,  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in NP4", angle.x = 45, remove.isolate = T)
gg1 + gg2




# get tables ---------------------------------------------------------

rm(list=ls())

files <- list.files(path="cellchat.2/", pattern = "RData")
names <- str_sub(files, 10, -7)

for(i in 1:length(files)){
  load(paste0("cellchat.2/", files[i]))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  df.net <- subsetCommunication(cellchat)
  df.net.2 <- subsetCommunication(cellchat, slot.name = "netP") # signaling pathways
  write.csv(df.net, file = paste0("cellchat.2/net_", names[i] , ".csv"))
  write.csv(df.net.2, file = paste0("cellchat.2/net2_", names[i] , ".csv"))
  saveRDS(cellchat, file = paste0("cellchat.2/cellchat_", names[i], ".rds"))
}



# spatial imaging ---------------------------------------------------------

cellchat.1 <- readRDS("cellchat.2/cellchat_P034_Pre.rds")
cellchat.2 <- readRDS("cellchat.2/cellchat_P034_Post.rds")
cellchat.3 <- readRDS("cellchat.2/cellchat_P039_Pre.rds")
cellchat.4 <- readRDS("cellchat.2/cellchat_P039_Post.rds")

pathways.show <- c("MHC-I") 
#pathways.show <- c("APP") 
#pathways.show <- c("TNF") 

levels(cellchat.3@idents)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.4, 
                    signaling = pathways.show, 
                    layout = "spatial", 
#                    sources.use = 4,
#                    targets.use = 1,
                    idents.use = c(4),
                    edge.width.max = 2, 
                    alpha.image = 0.3, 
                    #                    vertex.weight = "incoming", 
                    vertex.size.max = 3, 
big.gap = 15, small.gap = 15,
                    vertex.label.cex = 5,
remove.isolate = F
)


netVisual_aggregate(cellchat.2, 
                    signaling = pathways.show, 
                    layout = "spatial",
#                    idents.use = c(4),
                    sources.use = 4,
targets.use = 1,
vertex.size.max = 3,
                    vertex.label.cex = 5)




# end ---------------------------------------------------------------------
sessionInfo()






