
# intro -------------------------------------------------------------------

# scRNAseq data, NP137 manuscript
# preprocessing


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd()

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
  library(scales)
  
})

set.seed(170323)



# loading from barcode matrices  -----------------------------------------------------------------

dirs <- list.dirs(path = "./data/scRNAseq/matrices", recursive = F, full.names = T)

names <- str_sub(dirs, -3, -1)

list.samples.0 <- list()

for(i in 1:length(dirs)){
  list.samples.0[[i]] <- Read10X(data.dir = dirs[i])
  list.samples.0[[i]] <- CreateSeuratObject(list.samples.0[[i]])
  list.samples.0[[i]]$sample <- names[[i]]
}

names(list.samples.0) <- names
lapply(list.samples.0, dim)



# doublet estimation -------------------------------------------------------

object.doublet <- list.samples.0
for(i in 1:length(object.doublet)){
  
  temp1 <- SCTransform(object.doublet[[i]])
  temp1 <- RunPCA(temp1)
  temp1 <- RunUMAP(temp1, dims = 1:10)
  
  nExp_poi <- round(0.075*nrow(temp1@meta.data))  ## Assuming 7.5% doublet formation rate
  
  object.doublet[[i]] <- doubletFinder_v3(temp1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  rm(temp1)
  
}

#lapply(object.doublet, function(x) table(x@meta.data[, 16])) # for loom files
lapply(object.doublet, function(x) table(x@meta.data[, 8]))
#$NP3
#Doublet Singlet 
#929   11456 
#$NP4
#Doublet Singlet 
#763    9405 
doublet.table <- as.data.frame(lapply(object.doublet, function(x) table(x@meta.data[, 8])))
write.csv(doublet.table, file = "doublet.table.csv")

doublet.plot <- function(x){
  DimPlot(x, group.by = colnames(x@meta.data)[8]) + ggtitle(x@meta.data$sample)
}

plist <- lapply(object.doublet, doublet.plot)
do.call(grid.arrange, plist)


# to be checked retrospectively:

save(object.doublet, file = "data/objects/object.doublet.RData")



# merging -----------------------------------------------------------------

rm(list=ls())

load("data/objects/object.doublet.RData")

for(i in 1:length(object.doublet)){
  temp1 <- object.doublet[[i]]
  colnames(temp1@meta.data)[8] = "doublets"
  object.doublet[[i]] <- temp1
}

sc.merged <- merge(object.doublet[[1]], y = object.doublet[[2]], 
                   add.cell.ids = c("NP3","NP4"), project = "NP137")

table(sc.merged$doublets, sc.merged$sample)

sc.merged <- PercentageFeatureSet(sc.merged, pattern = "^MT-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.ribo")

VlnPlot(sc.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(sc.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        group.by = "doublets", ncol = 4)
plot(sc.merged$percent.mt, sc.merged$percent.ribo)

sc.merged <- subset(sc.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & doublets == "Singlet" & percent.mt < 25) %>%
  SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30) %>% RunTSNE()
# 16375 cells left

head(sc.merged[[]])
table(sc.merged$sample)
# NP3  NP4 
# 9216 7159

levels(x = sc.merged)

VlnPlot(sc.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
DimPlot(sc.merged, label = TRUE, label.size = 6)#, pt.size = 1.5)#, cols = "glasbey")

p1 <- DimPlot(sc.merged, reduction = "umap", group.by = "sample")#, pt.size = 1)
p2 <- DimPlot(sc.merged, reduction = "umap", label = F)#, cols = "glasbey")#, pt.size = 1)
plot_grid(p1, p2)

DimPlot(sc.merged, reduction = "umap", split.by = "sample", label = T, label.size = 4)#, cols = "glasbey")#, pt.size = 1)

# proportions

table(paste0(sc.merged@meta.data$sample, Idents(sc.merged)))
tab1 <- table(sc.merged@meta.data$sample, Idents(sc.merged))

write.csv(tab1, "results/scRNAseq/seurat/all.cells/initial_cluster_proportions_by_condition.csv")

save(sc.merged, file="data/objects/sc.merged.RData")



# cell cycle adjustment ---------------------------------------------------

rm(list=ls())

load(file="data/objects/sc.all.celltypes.RData")

DefaultAssay(sc.merged) <- "SCT"

RunPCA(sc.merged, features = VariableFeatures(sc.merged), ndims.print = 1:10, nfeatures.print = 10)
# no cell cycle genes in the first 10 PCs

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sc.merged <- CellCycleScoring(sc.merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(sc.merged[[]])
RidgePlot(sc.merged, features = c("MKI67", "TOP2A", "BIRC5", "MCM6"), ncol = 4)

# regress out cell cycle scores
sc.merged <- ScaleData(sc.merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sc.merged))

sc.merged <- RunPCA(sc.merged, 
                    features = VariableFeatures(sc.merged), 
                    ndims.print = 1:10, nfeatures.print = 10) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5)

pdf("results/scRNAseq/seurat/all.cells/Dimplot.minor.cell.types.cellcycle.pdf", width = 7, height = 6)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot.minor.cell.types.cellcycle.jpeg", 
     width = 800, height = 600)
DimPlot(sc.merged, label = T, group.by = "cell_type")
dev.off()


# no significant change in clusters after regressing cell cycle scores

save(sc.merged, file = "data/objects/sc.all.cellcycle.regressed.RData")




# end ---------------------------------------------------------------------
sessionInfo()






