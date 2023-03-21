
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



# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")

suppressPackageStartupMessages({
  library(Seurat);#  library(sctransform);  library(SeuratDisk); library(SeuratWrappers); library(SingleCellExperiment)
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



# inspection --------------------------------------------------------------

list.files("Nicolas/RDS")

immune <- readRDS("Nicolas/RDS/obj_patient_immune_sub_cluster_2_compute.RDS")
head(immune[[]])
table(immune$orig.ident)
tumoral <- readRDS("Nicolas/RDS/obj_patient_tumoral_sub_cluster_2_computed.RDS")
head(tumoral[[]])
table(tumoral$orig.ident)


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


object.list <- list.samples.0

object.list <- lapply(X = object.list, FUN = PercentageFeatureSet, pattern = "^mt-", col.name = "percent.mt") # not working for loom-based seurat object
object.list <- lapply(X = object.list, FUN = PercentageFeatureSet, pattern = "^Rp[sl][[:digit:]]", col.name = "percent.ribo")
#object.list <- lapply(X = object.list, FUN = subset, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)




# doublet estimation -------------------------------------------------------

object.doublet <- object.list
for(i in 1:length(object.doublet)){
  
  temp1 <- SCTransform(object.doublet[[i]])
  temp1 <- RunPCA(temp1)
  temp1 <- RunUMAP(temp1, dims = 1:10)
  
  nExp_poi <- round(0.075*nrow(temp1@meta.data))  ## Assuming 7.5% doublet formation rate
  
  object.doublet[[i]] <- doubletFinder_v3(temp1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  rm(temp1)
  
}

#lapply(object.doublet, function(x) table(x@meta.data[, 16])) # for loom files
lapply(object.doublet, function(x) table(x@meta.data[, 10]))
#$NP3
#Doublet Singlet 
#929   11456 
#$NP4
#Doublet Singlet 
#763    9405 
doublet.table <- as.data.frame(lapply(object.doublet, function(x) table(x@meta.data[, 10])))
write.csv(doublet.table, file = "doublet.table.csv")

doublet.plot <- function(x){
  DimPlot(x, group.by = colnames(x@meta.data)[10]) + ggtitle(x@meta.data$sample)
}

plist <- lapply(object.doublet, doublet.plot)
do.call(grid.arrange, plist)


# to be checked retrospectively:

save(object.doublet, file = "object.doublet.RData")



# merging -----------------------------------------------------------------


sc.merged <- merge(object.doublet[[1]], y = object.doublet[[2]], 
                   add.cell.ids = c("NP3","NP4"), project = "NP137")

sc.merged <- PercentageFeatureSet(sc.merged, pattern = "^MT-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.ribo") %>%
  subset(subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 25) %>%
  SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30) %>% RunTSNE()
# 11445 cells left

head(sc.merged[[]])
table(sc.merged$sample)
# NP3  NP4 
# 6111 5334

levels(x = sc.merged)

VlnPlot(sc.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
DimPlot(sc.merged, label = TRUE)#, pt.size = 1.5)#, cols = "glasbey")

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

known.markers <- c("PTPRC","PECAM1","EPCAM","PGR","TFF3","ACTA2")
FeaturePlot(sc.merged, known.markers, ncol = 3)
FeaturePlot(sc.merged, "PTPRC", ncol = 1)
VlnPlot(sc.merged, "CD19", ncol = 1) 
# 1 3 4 8 9 13 15 17 18 20 22
plot_density(sc.merged, reduction = "umap", known.markers)


## SciBet ----

DefaultAssay(sc.merged)
query <- GetAssayData(sc.merged, slot = "data")
query <- as.matrix(query)
query <- t(query)
query[1:10, 1:10]

# download models here: http://scibet.cancer-pku.cn/download_references.html
#model <- readr::read_csv("~/Dropbox/HH/code/references/scibet_major_human_cell_types.csv")
model <- readr::read_csv("~/Dropbox/HH/code/references/GSE115978_normal_scibet_core.csv")
#model <- readr::read_csv("~/Dropbox/HH/code/references/GSE115978_cancer_scibet_core.csv")
# Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq:
model <- readr::read_csv("~/Dropbox/HH/code/references/GSE72056_scibet_core.csv")

model <- pro.core(model)
colnames(model)
prd <- LoadModel(model) # generate Bet function from a model matrix

label.R <- prd(query)
head(label.R)
table(label.R)

sc.merged <- AddMetaData(sc.merged, label.R, "Scibet_1")
sc.merged <- AddMetaData(sc.merged, label.R, "Scibet_2")

head(sc.merged[[]])

DimPlot(sc.merged, group.by = "Scibet_1", label = T, repel = T) + NoLegend()
DimPlot(sc.merged, group.by = "Scibet_2", label = T, repel = T) + NoLegend()


## SingleR ----

ref <- HumanPrimaryCellAtlasData()
ref

table(ref$label.main)
length(table(ref$label.main)) # 36 cell types
table(ref$label.fine)
length(table(ref$label.fine)) # 157 cell types

pred.R <- SingleR(test = as.SingleCellExperiment(sc.merged), ref = ref, labels = ref$label.main)
pred.R <- SingleR(test = as.SingleCellExperiment(sc.merged), ref = ref, labels = ref$label.fine)
pred.R
table(pred.R$labels)

par(mfrow = c(1,1), mar=c(10,4,4,4))
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "Normal cell subtypes \n(SingleR - main)", border = 0)
par(mfrow = c(1,1), mar=c(22,4,4,4))
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "Normal cell subtypes \n(SingleR - fine)", border = 0)

sc.merged <- AddMetaData(sc.merged, pred.R$labels, "SingleR_main")
sc.merged <- AddMetaData(sc.merged, pred.R$labels, "SingleR_fine")
head(sc.merged[[]])

DimPlot(sc.merged, group.by = "SingleR_main", label = T, repel = T) + NoLegend()
DimPlot(sc.merged, group.by = "SingleR_fine", label = T, repel = T) + NoLegend()

Idents(sc.merged)

sc.merged <- RenameIdents(object = sc.merged, 
                          '0' = 'Epithelial', 
                          '5' = 'Epithelial', 
                          '12' = 'Epithelial', 
                          '16' = 'Epithelial', 
                          '21' = 'Epithelial', 
                          '23' = 'Epithelial', 
                          '7' = 'Endothelial', 
                          '11' = 'CAFs', 
                          '15' = 'Immune', 
                          '1' = 'Immune', 
                          '17' = 'Immune', 
                          '4' = 'Immune', 
                          '9' = 'Immune', 
                          '10' = 'Immune', 
                          '6' = 'Immune', 
                          '14' = 'Immune', 
                          '2' = 'Immune', 
                          '20' = 'Immune', 
                          '22' = 'Immune', 
                          '18' = 'Immune', 
                          '19' = 'Immune', 
                          '13' = 'Immune', 
                          '8' = 'Immune', 
                          '3' = 'Immune' 
                          )

table(Idents(sc.merged), sc.merged$sample)
#NP3  NP4
#Epithelial  2719  661
#Endothelial  373  152
#CAFs         233   73
#Immune      2786 4448

save(sc.merged, file="sc.all.celltypes.RData")

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
DimPlot(immune, group.by = "Scibet_2")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_main")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_fine")#, pt.size = 1)


# ScType
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


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
  ggtitle("Predicted Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2

DimPlot(immune, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'customclassif', cols = cluster.colors)        
DimPlot(immune, reduction = "umap", label = F, 
        repel = TRUE, group.by = 'customclassif',
        cols = cluster.colors,
        split.by = "sample")


save(immune, file = "immune.RData")




# T cells -------------------------------------------------

rm(list=ls())

load("immune.RData")

Idents(immune) <- immune$customclassif
table(immune$customclassif)

Tcells <- subset(immune, idents = c('Memory CD8+ T cells',
                                    'Naive CD4+ T cells',
                                    'Natural killer  cells',
                                    'CD8+ NKT-like cells'
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



# ScType
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# get cell-type-specific gene sets from our in-built database (DB)
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare(db_, tissue)
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max = sctype_score(scRNAseqData = Tcells[["SCT"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(Tcells@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Tcells@meta.data[Tcells@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Tcells@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


# plot

Tcells@meta.data$customclassif2 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Tcells@meta.data$customclassif2[Tcells@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(Tcells, reduction = "umap", label = TRUE, label.size = 3, label.box = T, repel = TRUE, group.by = 'customclassif2')
DimPlot(Tcells, reduction = "umap", label = F, repel = TRUE, group.by = 'customclassif2', split.by = "sample")



## ProjecTILs ----

# https://carmonalab.github.io/ProjecTILs_CaseStudies/
# https://github.com/carmonalab/ProjecTILs
# https://www.nature.com/articles/s41467-021-23324-4

load("Tcells.RData")

ref <- load.reference.map() # corresponds to mouse TILs
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



## new categories ----

rm(list=ls())

load("Tcells.RData")

Tcells <- querydata
rm(querydata)
head(Tcells[[]])

# automatic classifications

DimPlot(Tcells, group.by = "customclassif") # sctype
DimPlot(Tcells, group.by = "customclassif2") # sctype
DimPlot(Tcells, group.by = "functional.cluster") # ProjectTILs
DimPlot(Tcells, group.by = "Scibet_1") # 
DimPlot(Tcells, group.by = "Scibet_2") # 
DimPlot(Tcells, group.by = "SingleR_main") # 
DimPlot(Tcells, group.by = "SingleR_fine", 
        repel = T, label = T) + NoLegend() # 

# known markers

cat.colors <- as.vector(alphabet2())
DimPlot(Tcells, cols = cat.colors[1:10])

FeaturePlot(Tcells, "PTPRC")
FeaturePlot(Tcells, "CD3E")
FeaturePlot(Tcells, "CD4")
FeaturePlot(Tcells, "CD8A")

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
#    "KLRD1","NKG7","GNLY","CD3E"
  "NCAM1","CD14","CD19","CD3E"
)) # NK

# exhaustion: https://www.nature.com/articles/s41467-022-35238-w
plot_density(Tcells, reduction = "umap", c(
#  "CXCL13","CD274","PDCD1","PDCD1LG2"
  "CTLA4","CD80","CD86","CD28"
))


#library(UCell)
#signaturesHumanCellTypes <- readRDS("aux/signaturesHumanCellTypes.rds")

# CD8_Tpex
plot_density(Tcells, reduction = "umap", c(
#  "CD8A","CD8B",
#  "LAG3","XCL1","CRTAM", "TOX"
  "ZEB2",  "PDCD1", "TCF7", "CCR7" 
))

# CD8_Tex
plot_density(Tcells, reduction = "umap", c(
#  "CD8A", "CD8B", "LAG3", "PDCD1"
   "HAVCR2", "GZMB", "PRF1", "TIGIT" 
))


#CD3E: only in T Cells
#FCGR3A (CD16): in CD16+ monocytes and some expression NK cells
#GNLY: NK cells
#MS4A1: B cells
#GZMH: in GZMH+ T8 cells and some expression NK cells
#CD8A: in T8 cells
#CD4: in T4 and some myeloid cells
#CCR7: expressed more in memory cells 
#CD14: in CD14+ monocytes
#CD68: in monocytes/MF
#
#


Tcells <- RenameIdents(object = Tcells, 
                          '0' = 'memory CD4', 
                          '1' = 'CD8 effector memory', 
                          '2' = 'CD4', # resting T cell, circulating memory T cell
                          '3' = 'Naive CD4', 
                          '4' = 'CD4', # dividing T cell
                          '5' = 'Tregs', 
                          '6' = 'NK', # gamma delta ?
                          '7' = 'NK', # gamma delta ?
                          '8' = 'CD8 effector memory',
                          '9' = 'CD4'
)


head(Tcells[])
DimPlot(Tcells, cols = cat.colors[3:13])
DimPlot(Tcells, cols = cat.colors[3:13], split.by = "sample")


save(Tcells, file = "Tcells.RData")


#

rm(list=ls())

load("Tcells.RData")

Tcells$major.cats <- Idents(Tcells)
Idents(Tcells) <- Tcells$SCT_snn_res.0.8
DimPlot(Tcells, label = T, label.box = T)

load(file="~/Dropbox/HH/code/references/curated.ref.list.markers.Tcells.human.RData")

Tcells <- AddModuleScore_UCell(Tcells, features = ref.list)
signature.names <- paste0(names(ref.list), "_UCell")

VlnPlot(Tcells, features = signature.names[188:193])


Tcells <- RenameIdents(object = Tcells, 
                       '0' = 'memory CD4', 
                       '1' = 'CD8 effector memory', 
                       '2' = 'CD4 resting T cell', # resting T cell, circulating memory T cell
                       '3' = 'Naive CD4', 
                       '4' = 'CD4 dividing T cell', # dividing T cell
                       '5' = 'Tregs', 
                       '6' = 'NK like', # gamma delta ?
                       '7' = 'NK like', # gamma delta ?
                       '8' = 'CD8 effector memory',
                       '9' = 'Naive CD4'
)

cat.colors <- as.vector(alphabet2())

head(Tcells[])
library(scales)
show_col(cat.colors)
DimPlot(Tcells, cols = cat.colors[2:13])
DimPlot(Tcells, cols = cat.colors[2:13], split.by = "sample")


save(Tcells, file = "Tcells.RData")



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


# epithelial preprocessing ----

rm(list=ls())

load("epithelial.RData")

epithelial <- SCTransform(epithelial, vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)

p1 <- DimPlot(epithelial, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(epithelial, group.by = "sample")
p1+p2

DimPlot(epithelial, split.by = "sample")

# from: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

# related signatures from MsigDb:
library(msigdbr)
gsets_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gsets_msigdb <- gsets_msigdb[gsets_msigdb$gs_name==
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt.genes <- gsets_msigdb$gene_symbol

epithelial <- AddModuleScore(epithelial, 
                             name = "EMT_score",
                             list(emt.genes))

# compare to UCell:
library(UCell)
emt.genes <- list(EMT = emt.genes)
epithelial <- AddModuleScore_UCell(epithelial, features = emt.genes)
signature.names <- paste0(names(emt.genes), "_UCell")

VlnPlot(epithelial, features = "signature.names", group.by = "sample")
VlnPlot(epithelial, features = "EMT_score1", group.by = "sample")
VlnPlot(epithelial, 
        features = c("EMT_score1","EMT_UCell"),
        group.by = "sample")

head(epithelial[[]])
plot_density(epithelial, "EMT_UCell", reduction = "umap")
plot_density(epithelial, "EMT_score1", reduction = "umap")

FeaturePlot(epithelial, "EMT_UCell", split.by = "sample")
FeaturePlot(epithelial, "EMT_score1", split.by = "sample")

save(epithelial, file = "epithelial.RData")



# CellChat  ------------------------------------------

## grouped ----

rm(list=ls())

load("epithelial.RData")
load("Tcells.RData")

DimPlot(epithelial)
DimPlot(Tcells)

Idents(epithelial) <- "epithelial"
epithelial$cellchat.cat <- Idents(epithelial)
Idents(Tcells)
Tcells$cellchat.cat <- Idents(Tcells)
DefaultAssay(Tcells) <- "SCT"

sc.merged <- merge(epithelial, Tcells)
levels(sc.merged)
head(sc.merged[[]])

library(CellChat)
cellchat <- createCellChat(object = sc.merged, 
                           group.by = "cellchat.cat")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 8)
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

save(cellchat, file = "cellchat.RData")
save(sc.merged, file="sc.merged.RData")


## comparative ----

rm(list=ls())

load("sc.merged.RData")

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

jpeg("NP4_netVisual_circle_all.cells.jpeg", height = 900, width = 900, quality = 100)
netVisual_circle(cellchat.2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2, slot.name = "netP")
jpeg("NP4_netAnalysis_heatmap_all.cells.jpeg", height = 1161, width = 1065, quality = 100)
netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "outgoing", font.size = 10, height = 25) +
  netAnalysis_signalingRole_heatmap(cellchat.2, pattern = "incoming", font.size = 10, height = 25)
dev.off()

save(cellchat.2, file = "cellchat.NP4.RData")


#

rm(list=ls())

load("cellchat.NP3.RData")
load("cellchat.NP4.RData")

object.list <- list('NP3' = cellchat.1, 'NP4' = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

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
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "epithelial", comparison = c(1, 2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8 effector memory", comparison = c(1, 2))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK", comparison = c(1, 2))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3))


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
                                        title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(2, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
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

gg1 <- netVisual_bubble(cellchat, sources.use = c(1,2,4:7), targets.use = 3,  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in NP4", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2,4:7), targets.use = 3,  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in NP4", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4:7),  
                        comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in NP4", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4:7),  
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
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1,2,4:7), font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1,2,4:7), font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,2,4:7), targets.use = 3, font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1,2,4:7), targets.use = 3, font.size = 12,
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1,2,4:7), targets.use = 3, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 3, targets.use = c(1,2,4:7), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


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







