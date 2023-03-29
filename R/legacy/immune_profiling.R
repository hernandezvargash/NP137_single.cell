
# intro -------------------------------------------------------------------

# immune profiling with scRNAseq data
# revision to NP137 manuscript

# Referees' comment:
# There is a significant change in the immune landscape 
# in tumors treated with NP137 and 
# the authors also mentioned that EMT status is an important factor 
# contributing to resistance to immune-checkpoint inhibitors. 
# It is thus interesting to see whether NP137 can synergize 
# with current immune checkpoint therapies like anti-PD-1 or 
# even sensitize non-responding tumors to immune checkpoint blockade. 
# Before that, a more thorough analysis on their single cell data 
# is important, because they may discover altered pathways 
# in different immune cell population, and 
# this will be useful in guiding their choice of checkpoint inhibitors.

# tasks:
# identification of immune cell subtypes
# EMT score
# pathway analyses unbiased and targeted to immune checkpoints
# cell:cell communication

# 220323
# additional tasks to answer reviewer 3
# ECM score in tumor clusters and cellchat
# check T cell labelling (e.g. https://www.celltypist.org/)
# M1 vs M2
# check endothelial
# CAF subtypes
# checkpoint molecule expression

# also spatial analyses with cellchat



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
  
})

set.seed(170323)



# loading from barcode matrices  -----------------------------------------------------------------

dirs <- list.dirs(path = "./Nicolas/matrices", recursive = F, full.names = T)

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

save(object.doublet, file = "object.doublet.RData")



# merging -----------------------------------------------------------------

rm(list=ls())

load("object.doublet.RData")

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

write.csv(tab1, "initial_cluster_proportions_by_condition.csv")

save(sc.merged, file="sc.merged.RData")



# automated identification ------------------------------------------------

rm(list=ls())

load("sc.merged.RData")

head(sc.merged[[]])

known.markers <- c("PTPRC","PECAM1","EPCAM","PGR","TFF3","ACTA2","CD19","CD8A","CD4")
FeaturePlot(sc.merged, known.markers, ncol = 3)
plot_density(sc.merged, reduction = "umap", known.markers)


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
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
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

sc.merged@meta.data$ScType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc.merged@meta.data$ScType[sc.merged@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'ScType') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = F, repel = TRUE, group.by = 'ScType', split.by = "sample")


## relabel ----

Idents(sc.merged)
DimPlot(sc.merged, label = T, label.box = T)
DimPlot(sc.merged, label = T, split.by = "sample") + NoLegend()

sc.merged <- RenameIdents(object = sc.merged, 
                          '0' = 'Immune', 
                          '1' = 'Epithelial', 
                          '2' = 'Epithelial', 
                          '3' = 'Immune', 
                          '4' = 'Immune', 
                          '5' = 'Immune', 
                          '6' = 'Immune', 
                          '7' = 'Immune', 
                          '8' = 'Epithelial', 
                          '9' = 'Immune', 
                          '10' = 'Endothelial', 
                          '11' = 'Epithelial', 
                          '12' = 'Immune', 
                          '13' = 'CAFs', 
                          '14' = 'Immune', 
                          '15' = 'Immune', 
                          '16' = 'Immune', 
                          '17' = 'Immune', 
                          '18' = 'Epithelial', 
                          '19' = 'Immune', 
                          '20' = 'Immune', 
                          '21' = 'Endothelial', 
                          '22' = 'Immune', 
                          '23' = 'Immune'
                          )

sc.merged$cell_class <- Idents(sc.merged)

head(sc.merged[[]])
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'Scibet_1') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'Scibet_2') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'SingleR_main') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'SingleR_fine') + NoLegend()        
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'ScType') + NoLegend()        

Idents(sc.merged) <- sc.merged$seurat_clusters
sc.merged <- RenameIdents(object = sc.merged, 
                          '0' = 'CD4 T cells', 
                          '1' = 'Tumor cells', 
                          '2' = 'Tumor cells', 
                          '3' = 'CD8 NKT-like cells', # or CD8 T cells
                          '4' = 'Neutrophils', 
                          '5' = 'B cells', 
                          '6' = 'Macrophages', 
                          '7' = 'CD8 NKT-like cells', # or NK
                          '8' = 'Tumor cells', 
                          '9' = 'B cells', 
                          '10' = 'Endothelial', 
                          '11' = 'Tumor cells', 
                          '12' = 'B cells', 
                          '13' = 'CAFs', 
                          '14' = 'Basophils', 
                          '15' = 'B cells', 
                          '16' = 'Dendritic cells', 
                          '17' = 'Monocytes', 
                          '18' = 'Tumor cells', 
                          '19' = 'CD4 T cells', 
                          '20' = 'B cells', 
                          '21' = 'Tumor cells', 
                          '22' = 'Macrophages', 
                          '23' = 'B cells'
)

sc.merged$cell_type <- Idents(sc.merged)
head(sc.merged[[]])
DimPlot(sc.merged, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE) + NoLegend()        


table(Idents(sc.merged), sc.merged$sample)

save(sc.merged, file="sc.all.celltypes.RData")
#load(file="sc.all.celltypes.RData")


## proportions ----

tab0 <- table(sc.merged$sample, sc.merged$cell_type)
write.csv(tab0, file = "celltype.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:12]
condition.colors = c("brown","orange")

tab1 <- table(sc.merged$sample, sc.merged$seurat_clusters)
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cluster", "Sample","Proportion")
head(tab1)
p1 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  #  scale_fill_manual(values = cluster.colors) +
  ggtitle("Seurat Cluster Proportions") +
  theme_ipsum() +
  xlab("")

tab2 <- table(sc.merged$sample, sc.merged$cell_type)
tab2 <- as.data.frame(t(tab2))
colnames(tab2) <- c("Cell_Type", "Sample","Proportion")
head(tab2)
p2 <- ggplot(tab2, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

tab3 <- table(sc.merged$sample, sc.merged$cell_class)
tab3 <- as.data.frame(t(tab3))
colnames(tab3) <- c("Cell_Class", "Sample","Proportion")
head(tab3)
p3 <- ggplot(tab3, aes(x = Sample, y = Proportion, fill = Cell_Class)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = viridis(4)) +
    ggtitle("Cell Class Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2+p3


## create subsets ----

Idents(sc.merged) <- sc.merged$cell_class
immune <- subset(sc.merged, idents = "Immune")
epithelial <- subset(sc.merged, idents = "Epithelial")
endothelial <- subset(sc.merged, idents = "Endothelial")
CAFs <- subset(sc.merged, idents = "CAFs")

save(immune, file="immune.RData")
save(epithelial, file="epithelial.RData")
save(endothelial, file="endothelial.RData")
save(CAFs, file="CAFs.RData")



# immune cell identification ---------------------------------------------

rm(list=ls())

load(file="immune.RData")

immune <- SCTransform(immune, vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)

p1 <- DimPlot(immune, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(immune, group.by = "sample")
p1+p2

DimPlot(immune, split.by = "sample")#, pt.size = 1)
head(immune[[]])

DimPlot(immune, group.by = "Scibet_1")#, pt.size = 1)
DimPlot(immune, group.by = "Scibet_2")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_main")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_fine")#, pt.size = 1)
DimPlot(immune, group.by = "ScType")#, pt.size = 1)
DimPlot(immune, group.by = "cell_type")#, pt.size = 1)


# ScType
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
names(gs_list)

# get cell-type-specific gene sets from our in-built database (DB)
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare(db_, tissue)
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max = sctype_score(scRNAseqData = immune[["SCT"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(immune@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(immune@meta.data[immune@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(immune@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


# plot

immune@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  immune@meta.data$customclassif[immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(immune, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'customclassif')        
DimPlot(immune, reduction = "umap", label = F, repel = TRUE, group.by = 'customclassif', split.by = "sample")


#plot_density(immune, "CD8A", reduction = "umap")

# remove potentially contaminating cancer cells
 
immune <- subset(immune, 
       idents = c('13','19'),
       invert = T)


# proportions

tab1 <- table(immune$sample, immune$seurat_clusters)
write.csv(tab1, file = "immune.seurat.recluster.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:12]
condition.colors = c("brown","orange")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cell_type", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
#  scale_fill_manual(values = cluster.colors) +
  ggtitle("Seurat Cluster Proportions") +
  theme_ipsum() +
  xlab("")


tab2 <- table(immune$sample, immune$customclassif)
write.csv(tab2, file = "immune.predicted.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("Cell_type", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Predicted Immune Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2

DimPlot(immune, reduction = "umap", 
        label = TRUE, label.size = 3, label.box = F, 
        repel = TRUE, group.by = 'customclassif', 
        cols = cluster.colors) +
  ggtitle("Immune Cell Subtypes")
DimPlot(immune, reduction = "umap", label = F, 
        repel = TRUE, group.by = 'customclassif',
        cols = cluster.colors,
        split.by = "sample") +
  ggtitle("Immune Cell Subtypes")


save(immune, file = "immune.RData")


# save for CellTypist:

DefaultAssay(immune) <- "RNA"
immune.log <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
SaveH5Seurat(immune.log, filename = "immune.h5Seurat")
Convert("immune.h5Seurat", dest = "h5ad")


# opening the adata object from celltypist is not working
# so just updating the metadata

meta.immune <- read.csv("metadata.immune.celltypist.csv", row.names = 1)
head(meta.immune)

# this is another option
library(zellkonverter)
sce1=readH5AD("immune_celltypist.h5ad", verbose = TRUE)
adata_Seurat <- as.Seurat(sce1, counts = "X", data = NULL)
head(adata_Seurat[[]])
DimPlot(adata_Seurat, label = T, repel = T,
        group.by = "majority_voting") + NoLegend()

# this didn't work

Convert("immune_celltypist.h5ad", dest = "h5seurat", overwrite = TRUE)
test <- LoadH5Seurat("immune_celltypist.h5seurat")
head(test[[]])

library("Seurat")
library("anndata")
print("Convert from Scanpy to Seurat...")
data <- read_h5ad("immune_celltypist.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
print(str(data))



# T cells -------------------------------------------------

rm(list=ls())

load("immune.RData")

Idents(immune) <- immune$customclassif
table(immune$customclassif)

Tcells <- subset(immune, idents = c('Memory CD8+ T cells',
                                    'Naive CD4+ T cells',
                                    'Natural killer  cells',
                                    'Naive CD8+ T cells'
                                    ), invert = F)
table(Tcells$sample, Tcells$customclassif)

# increased resolution
Tcells <- SCTransform(Tcells, vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% RunUMAP(dims = 1:30)

p1 <- DimPlot(Tcells, label = F, 
              group.by = "customclassif",
              label.box = F, label.size = 5)
p2 <- DimPlot(Tcells, group.by = "sample")
p1+p2

DimPlot(Tcells, 
        group.by = "customclassif",
        split.by = "sample")#, pt.size = 1)


save(Tcells, file = "Tcells.RData")



## ProjecTILs ----

# https://carmonalab.github.io/ProjecTILs_CaseStudies/
# https://github.com/carmonalab/ProjecTILs
# https://www.nature.com/articles/s41467-021-23324-4

load("Tcells.RData")

ref <- load.reference.map() # corresponds to mouse TILs, as no human available
head(ref[[]])
table(ref$functional.cluster)
table(ref$TILPRED)

DefaultAssay(Tcells) <- "RNA"

# projection:

#projected <- Run.ProjecTILs(sc.sub, ref=ref)
#head(projected[[]])
#plot.projection(ref, projected, linesize = 0.5, pointsize = 0.5)
#plot.statepred.composition(ref, projected, metric = "Percent")
#plot.states.radar(ref, query = projected, min.cells = 30)

# classification:

querydata <- ProjecTILs.classifier(query = Tcells, 
                                   ref = ref, 
                                   filter.cells = F,
                                   split.by = "sample")
head(querydata[[]])
table(querydata$functional.cluster)

p1 <- DimPlot(querydata, label = T, label.box = F, label.size = 3)
p2 <- DimPlot(querydata, group.by = "customclassif", label = T) + NoLegend()
p3 <- DimPlot(querydata, group.by = "functional.cluster", label = T) + NoLegend()

p2+p3

DimPlot(querydata, group.by = "functional.cluster", label = F, split.by = "sample")

head(querydata[[]])
save(querydata, file = "Tcells.RData")

# ProjecTILs wrongly re-labels NK cells,
# likely because this cell type is not present in the reference
# However, it identifies Tregs and Tex cells


## new categories ----

rm(list=ls())

load("Tcells.RData")

Tcells <- querydata
rm(querydata)
head(Tcells[[]])

# automatic classifications

DimPlot(Tcells, label = T)
DimPlot(Tcells, group.by = "customclassif") # sctype
DimPlot(Tcells, group.by = "functional.cluster", 
        cols = as.vector(glasbey(n=10)),
        label = T) # ProjectTILs
DimPlot(Tcells, group.by = "Scibet_1") # 
DimPlot(Tcells, group.by = "Scibet_2") # 
DimPlot(Tcells, group.by = "SingleR_main") # 
DimPlot(Tcells, group.by = "SingleR_fine", 
        repel = T, label = T) + NoLegend() # 

# top cluster (7) is identified as Tex by ProjecTILs and 
# as gammadelta by SingleR_fine


# known markers

cat.colors <- as.vector(alphabet2())
DimPlot(Tcells, cols = cat.colors[1:10])

plot_density(Tcells, 
             reduction = "umap",
             c("PTPRC","CD3E","CD4","CD8A"))
plot_density(Tcells, reduction = "umap",
             c("CCR7", "SELL", "IL7R", "CD27", "LEF1")) # CD8 naive
plot_density(Tcells, reduction = "umap", c(
#               "IL2RA", "TNFRSF8", "CD69", "TNFRSF4" 
#               "ICOS", "KLRG1", "HAVCR2", "TBX21"
#               "IFNG", "IL2", "PRF1", "GZMB"#, "GZMA"
               "TNF", "CCL3", "CCL4", "CCL5"
               )) # CD8 effector
plot_density(Tcells, reduction = "umap", c(
#               "CD44", "SELL", "IL7R", "KLRG1"
               "CCR7", "EOMES", "TBX21", "IFNG"
#               "TNF", "IL2", "PRF1","CD8A"
)) # CD8 effector memory: (EMOS, GZMK, IFNG)?
#  CD62L low, CD127 hi, KLRG1 hi, 
#  CCR7 lo, IFNG++, TNFA++, 
plot_density(Tcells, reduction = "umap", c(
  "KLRB1", "SLC4A1", "ZBTB16","CD8A"
)) # MAIT cells
plot_density(Tcells, reduction = "umap", c(
#  "CD28", "KLRG1","PDCD1", "IKZF2"
#  "LAG3", "HLA-DRB1", "B3GAT1","FOXP3"
  "IKZF1", "EGR1", "EGR2","CD8A"
)) # regulatory: CD28-, KLRG1++

plot_density(Tcells, reduction = "umap", c(
#  "CCR1", "CCR5", "STAT1","STAT4"
#  "KLRD1", "IFNGR1", "CXCR3","CXCR6"
  "IL2", "TNF","LTA", "IFNG" #"TBX21"
)) # Th1
plot_density(Tcells, reduction = "umap", c(
  "IL12RB1", "IL18R1", "TNFSF11","HAVCR2"
)) # effector
plot_density(Tcells, reduction = "umap", c(
  "CXCR4", "CCR4", "CCR8", "CD4" #,"PTGDR2"
#  "HAVCR1", "IL17RB", "IL33","BATF"
#  "GATA3", "IRF4", "STAT6","IL4"
#   "IL13", "AREG", "CCR4", "CXCR4" # "IL5"
)) # Th2
plot_density(Tcells, reduction = "umap", c(
#  "IL1R1", "KLRB1", "CCR4", "IL21R"
#  "IL12RB1", "AHR", "BATF","MAF"
#  "NFKBIZ", "IRF4", "RORA","RORC"
  "STAT3", "IL17F", "IL21", "IL22"
)) # Th17 # (CD4, RORC, IL23R, CSF2)?
plot_density(Tcells, reduction = "umap", c(
  "TBX21", "GZMB", "GZMH", "PRF1"
)) # cytotoxic
plot_density(Tcells, reduction = "umap", c(
  "IL2RA", "CTLA4", "FOXP3", "CD4"
)) # Tregs
plot_density(Tcells, reduction = "umap", c(
#  "CXCR3", "CXCR5", "ICOS", "PDCD1"
#  "BATF", "BCL6", "MAF", "IRF4"
#  "STAT3", "IL21", "CD4","ICOS"
    "CD40LG", "TOX2", "CD200",  "BATF"  
)) # TFH
rownames(Tcells)[grep("TRG", rownames(Tcells))]
rownames(Tcells)[grep("TRD", rownames(Tcells))]

FeaturePlot(Tcells, "TRGC1")
plot_density(Tcells, reduction = "umap", c(
  "TRGC1","TRGC2","TRGV8","TRDC"
)) # gamma-delta , https://www.pnas.org/doi/10.1073/pnas.1818488116
plot_density(Tcells, reduction = "umap", c(
    "KLRD1","NKG7","GNLY","CD3E"
#  "NCAM1","CD14","CD19","CD3E"
)) # NK

# exhaustion: https://www.nature.com/articles/s41467-022-35238-w
plot_density(Tcells, reduction = "umap", c(
  "CXCL13","CD274","PDCD1","PDCD1LG2"
#  "CTLA4","CD80","CD86","CD28"
))

# CD8_Tpex
plot_density(Tcells, reduction = "umap", c(
#  "CD8A","CD8B",
  "LAG3","XCL1","CRTAM", "TOX"
#  "ZEB2",  "PDCD1", "TCF7", "CCR7" 
))

# CD8_Tex
plot_density(Tcells, reduction = "umap", c(
#  "CD8A", "CD8B", "LAG3", "PDCD1"
   "HAVCR2", "GZMB", "PRF1", "TIGIT" 
))


Tex <- WhichCells(Tcells, 
                  expression = functional.cluster == "CD8_Tex" &
                    SCT_snn_res.0.8 == '7')
DimPlot(Tcells, cells.highlight = Tex)
#Tcells <- SetIdent(Tcells, cells = Tex, value = 'CD8 Tex')

gd <- WhichCells(Tcells, 
                  expression = SingleR_fine == "T_cell:gamma-delta" &
                    TRGC1 > 0)
DimPlot(Tcells, cells.highlight = gd)

Tcells <- RenameIdents(object = Tcells, 
                          '0' = 'CD4 T cells', 
                          '1' = 'CD8 T cells', 
                          '2' = 'CD4 T cells',
                          '3' = 'Natural killer cells', 
                          '4' = 'Treg cells',
                          '5' = 'Natural killer cells', 
                          '6' = 'CD4 T cells',
                          '7' = 'CD8 Tex',
                          '8' = 'Natural killer cells',
                          '9' = 'CD8 T cells',
                       '10' = 'Undefined',
                       '11' = 'CD4 T cells'
                       )


head(Tcells[])
DimPlot(Tcells, cols = cat.colors[3:13])
DimPlot(Tcells, cols = cat.colors[3:13], split.by = "sample")

Tcells$final.cats <- Idents(Tcells)

save(Tcells, file = "Tcells.RData")
meta <- Tcells@meta.data
rownames(meta) <- colnames(Tcells)
write.csv(meta, file="metadata.Tcells.csv")



## cell type proportions ---------------------------------------------------

load("Tcells.RData")

tab1 <- table(Tcells$sample, Idents(Tcells))
prop.table(tab1)
write.csv(tab1, file = "Tcell.proportions.csv")

cat.colors <- as.vector(alphabet2())
#colors <- as.vector(glasbey(n=24))
cluster.colors = cat.colors[3:13]
condition.colors = c("brown","orange")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cell_Type", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("")

p2 <- DimPlot(Tcells, cols = cat.colors[3:13], split.by = "sample")
p3 <- DimPlot(Tcells, cols = cat.colors[3:13])

p1+p3



# save for CellTypist:

Tcells
DefaultAssay(Tcells) <- "RNA"
Tcells.log <- NormalizeData(Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
SaveH5Seurat(Tcells.log, filename = "Tcells.h5Seurat")
Convert("Tcells.h5Seurat", dest = "h5ad")


# or run it in R using RIRA_classification package

Tcells <- RunCellTypist(
  seuratObj,
  modelName = "Immune_All_Low.pkl",
  pThreshold = 0.5,
  minProp = 0,
  useMajorityVoting = TRUE,
  mode = "prob_match",
  extraArgs = c("--mode", mode, "--p-thres", pThreshold, "--min-prop", minProp),
  assayName = Seurat::DefaultAssay(seuratObj),
  columnPrefix = NULL,
  maxAllowableClasses = 6,
  minFractionToInclude = 0.01,
  minCellsToRun = 200,
  maxBatchSize = 1e+05,
  retainProbabilityMatrix = FALSE,
  runCelltypistUpdate = TRUE
)

head(Tcells[[]])

# opening the adata object from celltypist is not working
# so just updating the metadata

meta.Tcells <- read.csv("metadata.Tcells.celltypist.csv", row.names = 1)
head(meta.Tcells)

# this is another option
library(zellkonverter)
sce1=readH5AD("Tcells_celltypist.h5ad", verbose = TRUE)
adata_Seurat <- as.Seurat(sce1, counts = "X", data = NULL)
head(adata_Seurat[[]])
DimPlot(adata_Seurat, label = T, repel = T,
        group.by = "majority_voting") + NoLegend()

# this didn't work

Convert("immune_celltypist.h5ad", dest = "h5seurat", overwrite = TRUE)
test <- LoadH5Seurat("immune_celltypist.h5seurat")
head(test[[]])

library("Seurat")
library("anndata")
print("Convert from Scanpy to Seurat...")
data <- read_h5ad("immune_celltypist.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
print(str(data))



# checkpoint genes --------------------------------------------------------

load("sc.all.celltypes.RData")
head(sc.merged[[]])
DimPlot(sc.merged, label = T)

known.markers <- c("CXCL9","CCR5","CXCL13",
                   "PDCD1","CD274",
                   "CTLA4","CD80","CD86")

jpeg("checkpoint.jpeg", width = 1200, height = 900, quality = 100)
FeaturePlot(sc.merged, known.markers, ncol = 3)
dev.off()

jpeg("checkpoint_density.jpeg", width = 1200, height = 900, quality = 100)
plot_density(sc.merged, reduction = "umap", known.markers)
dev.off()


i = 2
VlnPlot(sc.merged, known.markers[i], 
        group.by = "sample", ncol = 1) +
  stat_compare_means() + xlab("")
FeaturePlot(sc.merged, known.markers[i], 
            order = T, pt.size = 0.6,
            split.by = "sample")
VlnPlot(sc.merged, known.markers[i],
        group.by = "cell_type",
        split.by = "sample", ncol = 1)




# epithelial preprocessing ----

rm(list=ls())

load("epithelial.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("brown","orange")

epithelial <- SCTransform(epithelial, vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)

p1 <- DimPlot(epithelial, label = T, label.box = F, label.size = 5, cols = cluster.colors)
p2 <- DimPlot(epithelial, group.by = "sample", cols = cluster.colors)
p1+p2

DimPlot(epithelial, split.by = "sample", cols = cluster.colors)


prop.table(table(epithelial$sample, epithelial$seurat_clusters))

tab1 <- table(epithelial$sample, epithelial$seurat_clusters)
write.csv(tab1, file = "tumor.cluster.proportions.csv")
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cluster", "Sample","Proportion")
head(tab1)
ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Tumor Cluster Proportions") +
  theme_ipsum() +
  xlab("")

table(epithelial$sample)
# cluster 1 is 18% in NP3, compared to 13% in NP4 (see below)



# from: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

# related signatures from MsigDb:
library(msigdbr)
gsets_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gsets_msigdb <- gsets_msigdb[gsets_msigdb$gs_name==
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt.genes <- gsets_msigdb$gene_symbol


# score with AddModuleScore
epithelial <- AddModuleScore(epithelial, 
                             name = "EMT_score",
                             list(emt.genes))

v1 <- VlnPlot(epithelial, features = "EMT_score1", group.by = "sample") +
  ggtitle("EMT Score [Seurat]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))



# compare to UCell:
library(UCell)
emt.genes <- list(EMT = emt.genes)
epithelial <- AddModuleScore_UCell(epithelial, features = emt.genes)
signature.names <- paste0(names(emt.genes), "_UCell")

library(ggpubr)
v2 <- VlnPlot(epithelial, features = signature.names, group.by = "sample") +
  ggtitle("EMT Score [UCell]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))



# compare with AUCell
library(AUCell)
HF <- list(geneSet1=c("Cxcl14","Wfdc18","Sbsn","Shisa2","Krt17","Gja1","Postn","Selenop","Gjb2","Hspb1",
                      "Clasrp","Lhx2","Pthlh","Gas1","Sox9","Id4","Sgk1","Sdc1","Pdzrn3","Tbx1"))
cells_AUC_EMT <- AUCell_run(epithelial@assays$RNA@counts, emt.genes)
epithelial[["EMT_cells"]] <- cells_AUC_EMT@assays@data@listData[["AUC"]][1,]

v3 <- VlnPlot(epithelial,features = "EMT_cells", group.by = "sample") +
  ggtitle("EMT Score [AUCell]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))


grid.arrange(v1,v2,v3, ncol = 3)


VlnPlot(epithelial,features = "EMT_cells", 
        split.by = "sample",
        group.by = "seurat_clusters") +
  ggtitle("EMT Score [AUCell]") + 
  xlab("") +
#  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))



head(epithelial[[]])
d1 <- plot_density(epithelial, "EMT_UCell", reduction = "umap") +
  ggtitle("EMT Score [Seurat]")
d2 <- plot_density(epithelial, "EMT_score1", reduction = "umap") +
  ggtitle("EMT Score [UCell]")
d3 <- plot_density(epithelial, "EMT_cells", reduction = "umap") +
  ggtitle("EMT Score [AUCell]")

FeaturePlot(epithelial, 
            c("EMT_UCell","EMT_score1","EMT_cells"),
            reduction = "umap")

d1 + d2 + d3

DimPlot(epithelial, label = T, cols = cluster.colors, label.box = T)
DimPlot(epithelial, label = T, split.by = "sample")


save(epithelial, file = "epithelial.RData")




# end ---------------------------------------------------------------------







