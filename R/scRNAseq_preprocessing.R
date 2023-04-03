
# intro -------------------------------------------------------------------

# scRNAseq data, NP137 manuscript
# # preprocessing


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



# automated identification ------------------------------------------------

rm(list=ls())

load("data/objects/sc.merged.RData")

head(sc.merged[[]])

known.markers <- c("PTPRC","PECAM1","EPCAM","PGR","TFF3","ACTA2","CD19","CD8A","CD4")
pdf("results/scRNAseq/seurat/all.cells/FeaturePlot.known.markers.pdf", 
    width = 15, height = 12)
jpeg("results/scRNAseq/seurat/all.cells/FeaturePlot.known.markers.jpeg", 
    quality = 100, height = 900, width = 1200)
FeaturePlot(sc.merged, known.markers, ncol = 3)
dev.off()
pdf("results/scRNAseq/seurat/all.cells/DensityPlot.known.markers.pdf", 
    width = 15, height = 12)
jpeg("results/scRNAseq/seurat/all.cells/DensityPlot.known.markers.jpeg", 
     quality = 100, height = 900, width = 1200)
plot_density(sc.merged, reduction = "umap", known.markers)
dev.off()


## SciBet ----

DefaultAssay(sc.merged)
query <- GetAssayData(sc.merged, slot = "data")
query <- as.matrix(query)
query <- t(query)
query[1:10, 1:10]

# download models here: http://scibet.cancer-pku.cn/download_references.html
model <- readr::read_csv("~/Dropbox/HH/code/references/GSE115978_normal_scibet_core.csv")
model <- pro.core(model)
colnames(model)
prd <- LoadModel(model) # generate Bet function from a model matrix
label.R <- prd(query)
head(label.R)
table(label.R)
sc.merged <- AddMetaData(sc.merged, label.R, "Scibet_1")
head(sc.merged[[]])
rm(model)


model <- readr::read_csv("~/Dropbox/HH/code/references/GSE72056_scibet_core.csv")
model <- pro.core(model)
colnames(model)
prd <- LoadModel(model) # generate Bet function from a model matrix
label.R <- prd(query)
head(label.R)
table(label.R)
sc.merged <- AddMetaData(sc.merged, label.R, "Scibet_2")
head(sc.merged[[]])

DimPlot(sc.merged, group.by = "Scibet_1", label = T, label.size = 6, repel = T) + NoLegend()
DimPlot(sc.merged, group.by = "Scibet_2", label = T, label.size = 6, repel = T) + NoLegend()


## SingleR ----

ref <- HumanPrimaryCellAtlasData()
ref

table(ref$label.main)
length(table(ref$label.main)) # 36 cell types
table(ref$label.fine)
length(table(ref$label.fine)) # 157 cell types

pred.R <- SingleR(test = as.SingleCellExperiment(sc.merged), ref = ref, labels = ref$label.main)
par(mfrow = c(1,1), mar=c(10,4,4,4))
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "Normal cell subtypes \n(SingleR - main)", border = 0)
sc.merged <- AddMetaData(sc.merged, pred.R$labels, "SingleR_main")
rm(pred.R)

pred.R <- SingleR(test = as.SingleCellExperiment(sc.merged), ref = ref, labels = ref$label.fine)
par(mfrow = c(1,1), mar=c(22,4,4,4))
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "Normal cell subtypes \n(SingleR - fine)", border = 0)
sc.merged <- AddMetaData(sc.merged, pred.R$labels, "SingleR_fine")
head(sc.merged[[]])

DimPlot(sc.merged, group.by = "SingleR_main", label = T, repel = T) + NoLegend()
DimPlot(sc.merged, group.by = "SingleR_fine", label = T, repel = T) + NoLegend()


## ScType ----
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", "Lung") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max = sctype_score(scRNAseqData = sc.merged[["SCT"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(sc.merged@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sc.merged@meta.data[sc.merged@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc.merged@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


# plot

sc.merged@meta.data$ScType_Immune = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc.merged@meta.data$ScType_Immune[sc.merged@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'ScType_Immune') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = F, repel = TRUE, group.by = 'ScType_Immune', split.by = "sample")

sc.merged@meta.data$ScType_Lung = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc.merged@meta.data$ScType_Lung[sc.merged@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'ScType_Lung') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = F, repel = TRUE, group.by = 'ScType_Lung', split.by = "sample")

head(sc.merged[[]])
save(sc.merged, file="data/objects/sc.all.celltypes.RData")



# relabel ----

load(file="data/objects/sc.all.celltypes.RData")

# inspect different classifications


d1 <- DimPlot(sc.merged, group.by = "seurat_clusters", label = T, label.size = 6, repel = T, label.box = T) + NoLegend()
d1b <- DimPlot(sc.merged, group.by = "seurat_clusters", label = F, label.size = 6, repel = F, label.box = F)

d2 <- DimPlot(sc.merged, group.by = "Scibet_1", label = T, label.size = 6, repel = T) + NoLegend()
d3 <- DimPlot(sc.merged, group.by = "Scibet_2", label = T, label.size = 6, repel = T) + NoLegend()
d4 <- DimPlot(sc.merged, group.by = "SingleR_main", label = T, repel = T) + NoLegend()
d5 <- DimPlot(sc.merged, group.by = "SingleR_fine", label = T, repel = T) + NoLegend()
d6 <- DimPlot(sc.merged, group.by = "ScType_Immune", label = T, repel = T) + NoLegend()
d7 <- DimPlot(sc.merged, group.by = "ScType_Lung", label = T, repel = T) + NoLegend()

pdf("results/scRNAseq/seurat/all.cells/Dimplot.all.cells.pdf", width = 14, height = 7)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot.all.cells.jpeg", 
     width = 1200, height = 600, quality = 100)
d1 + d1b
dev.off()

pdf("results/scRNAseq/seurat/all.cells/Dimplots.predicted.cell.types.pdf", width = 18, height = 10)
jpeg("results/scRNAseq/seurat/all.cells/Dimplots.predicted.cell.types.jpeg", 
     width = 1300, height = 700, quality = 100)
grid.arrange(d2,d3,d4,d5,d6,d7, ncol = 3)
dev.off()

Idents(sc.merged) <- sc.merged$seurat_clusters

DimPlot(sc.merged, label = T, label.box = T)
DimPlot(sc.merged, label = T, split.by = "sample") + NoLegend()

sc.merged <- RenameIdents(object = sc.merged, 
                          '0' = 'Immune cells', 
                          '1' = 'Tumor cells', 
                          '2' = 'Tumor cells', 
                          '3' = 'Immune cells', 
                          '4' = 'Immune cells', 
                          '5' = 'Immune cells', 
                          '6' = 'Immune cells', 
                          '7' = 'Immune cells', 
                          '8' = 'Tumor cells', 
                          '9' = 'Immune cells', 
                          '10' = 'Endothelial cells', 
                          '11' = 'Tumor cells', 
                          '12' = 'Immune cells', 
                          '13' = 'CAFs', 
                          '14' = 'Immune cells', 
                          '15' = 'Immune cells', 
                          '16' = 'Immune cells', 
                          '17' = 'Immune cells', 
                          '18' = 'Tumor cells', 
                          '19' = 'Immune cells', 
                          '20' = 'Immune cells', 
                          '21' = 'Tumor cells', 
                          '22' = 'Immune cells', 
                          '23' = 'Immune cells'
)

sc.merged$cell_class <- Idents(sc.merged)

colors <- glasbey()
show_col(colors)
colors1 <- colors[c(8,15,7,21)]

pdf("results/scRNAseq/seurat/all.cells/Dimplot1.major.cell.types.pdf", width = 8, height = 6)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot1.major.cell.types.jpeg", 
     width = 800, height = 600)
DimPlot(sc.merged, cols = colors1)# + NoLegend()
dev.off()

Idents(sc.merged) <- sc.merged$sample
sc.merged <- RenameIdents(sc.merged, "NP3" = "C1D1",
                          "NP4" = "C3D1")
sc.merged$ID <- Idents(sc.merged)
Idents(sc.merged) <- sc.merged$cell_class

pdf("results/scRNAseq/seurat/all.cells/Dimplot2.major.cell.types.pdf", width = 12, height = 6)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot2.major.cell.types.jpeg", 
     width = 1200, height = 600)
DimPlot(sc.merged, cols = colors1, split.by = "ID")# + NoLegend()
dev.off()


# fine labeling

Idents(sc.merged) <- sc.merged$seurat_clusters
sc.merged <- RenameIdents(object = sc.merged, 
                          '0' = 'CD4+ T cells', 
                          '1' = 'Tumor cells', 
                          '2' = 'Tumor cells', 
                          '3' = 'CD8+ T cells', # or CD8 T cells
                          '4' = 'Neutrophils', 
                          '5' = 'B cells', 
                          '6' = 'Macrophages', 
                          '7' = 'NK cells', # or NK
                          '8' = 'Tumor cells', 
                          '9' = 'B cells', 
                          '10' = 'Endothelial cells', 
                          '11' = 'Tumor cells', 
                          '12' = 'B cells', 
                          '13' = 'CAFs', 
                          '14' = 'Basophils', 
                          '15' = 'B cells', 
                          '16' = 'Dendritic cells', 
                          '17' = 'Monocytes', 
                          '18' = 'Tumor cells', 
                          '19' = 'CD4+ T cells', 
                          '20' = 'B cells', 
                          '21' = 'Tumor cells', 
                          '22' = 'Macrophages', 
                          '23' = 'B cells'
)

sc.merged$cell_type <- Idents(sc.merged)
head(sc.merged[[]])

alphabet()
colors2 <- alphabet2()[1:12]
show_col(colors2)
names(colors2) <- levels(as.factor(sc.merged$cell_type))

d1 <- DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, cols = colors2) + NoLegend()        
d2 <- DimPlot(sc.merged, reduction = "umap", label = F, label.size = 3, label.box = T, repel = TRUE, cols = colors2) #+ NoLegend()        

pdf("results/scRNAseq/seurat/all.cells/Dimplot1.minor.cell.types.pdf", width = 12, height = 6)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot1.minor.cell.types.jpeg", 
     width = 1200, height = 600)
d1+d2
dev.off()

pdf("results/scRNAseq/seurat/all.cells/Dimplot2.minor.cell.types.pdf", width = 12, height = 6)
jpeg("results/scRNAseq/seurat/all.cells/Dimplot2.minor.cell.types.jpeg", 
     width = 1200, height = 500)
DimPlot(sc.merged, group.by = "cell_type", cols = colors2, split.by = "ID")# + NoLegend()
dev.off()


table(Idents(sc.merged), sc.merged$sample)

save(sc.merged, file="data/objects/sc.all.celltypes.RData")
#load(file="data/objects/sc.all.celltypes.RData")

meta <- sc.merged@meta.data
head(meta)
write.csv(meta, file = "data/objects/sc.all.metadata.csv")


# proportions ----

rm(list=ls())

meta <- read.csv(file = "data/objects/sc.all.metadata.csv")
head(meta)
meta$cell_class <- as.factor(meta$cell_class)
meta$cell_type <- as.factor(meta$cell_type)

colors <- glasbey()
colors1 <- colors[c(21,7,8,15)]
names(colors1) <- levels(as.factor(meta$cell_class))
colors2 <- alphabet2()[c(5,10,9,1,
                         3,11,8,6,
                         12,4,7,2)]
names(colors2) <- levels(as.factor(meta$cell_type))
colors3 = c("brown","orange")

tab1 <- table(meta$ID, meta$cell_class)
head(tab1)
write.csv(tab1, file = "results/scRNAseq/seurat/cell.class_proportions.csv")
tab2 <- table(meta$ID, meta$cell_type)
head(tab2)
write.csv(tab2, file = "results/scRNAseq/seurat/cell.type_proportions.csv")


tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cell.Class", "Sample","Proportion")
p1 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cell.Class)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors1) +
  ggtitle("Cell Class Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Cell Class Proportions") + NoLegend()
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cell.Class)) +
  geom_col() +
  scale_fill_manual(values = colors1) +
  ggtitle("Cell Class Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Cell Class Numbers")

png("results/scRNAseq/seurat/all.cells/Barplot.cellclass.png")
jpeg("results/scRNAseq/seurat/all.cells/Barplot.cellclass.jpeg", 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()


tab2 <- as.data.frame(t(tab2))
colnames(tab2) <- c("Cell.Type", "Sample","Proportion")
p1 <- ggplot(tab2, aes(x = Sample, y = Proportion, fill = Cell.Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors2) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Cell Type Proportions") + NoLegend()
p2 <- ggplot(tab2, aes(x = Sample, y = Proportion, fill = Cell.Type)) +
  geom_col() +
  scale_fill_manual(values = colors2) +
  ggtitle("Cell Type Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Cell Type Numbers")

png("results/scRNAseq/seurat/all.cells/Barplot.celltype.png")
jpeg("results/scRNAseq/seurat/all.cells/Barplot.celltype.jpeg", 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()


# create subsets ----

Idents(sc.merged) <- sc.merged$cell_class
immune <- subset(sc.merged, idents = "Immune cells") # 9885
tumor <- subset(sc.merged, idents = "Tumor cells") # 5622
endothelial <- subset(sc.merged, idents = "Endothelial cells") # 562
CAFs <- subset(sc.merged, idents = "CAFs") # 306

save(immune, file="data/objects/immune.RData")
save(tumor, file="data/objects/tumor.RData")
save(endothelial, file="data/objects/endothelial.RData")
save(CAFs, file="data/objects/CAFs.RData")



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






