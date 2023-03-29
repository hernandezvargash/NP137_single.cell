
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




# end ---------------------------------------------------------------------







