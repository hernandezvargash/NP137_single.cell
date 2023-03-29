
# immune profiling with spatial scRNAseq data
# revision to NP137 manuscript
# Visium avant (C1D1) et après traitement (C3D1)

# 220323, spatial tasks:
# ECM score in tumor clusters
# T cell identification
# identify other cell types: endothelial, CAF subtypes
# checkpoint molecule expression
# spatial analyses with cellchat

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
# https://satijalab.org/seurat/articles/spatial_vignette.html#working-with-multiple-slices-in-seurat-1
# http://giottosuite.com/index.html


# libraries ---------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(Seurat)
  #  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  
})


setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")



# preprocessing -----------------------------------------------------------


## loading ----

list.files()

samples <- list.dirs(path = "Nicolas/Visium/", full.names = F, recursive = F)

data.dir <- list.dirs(path = "Nicolas/Visium", full.names = T, recursive = F)

list.samples <- list()
for(i in 1:length(samples)){
  list.samples[[i]] <- Load10X_Spatial(data.dir = paste0(data.dir[i], "/outs/"))
  list.samples[[i]]$sample <- samples[[i]]
}

patient_ID <- str_sub(samples, 4, -6)
treatment <- c(rep(c("Pre","Post"), 2, by = 2))

names(list.samples) <- paste0("P", patient_ID, "_", treatment)
lapply(list.samples, dim)



## QC ----

list.samples <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^MT-", col.name = "percent.mt")
list.samples <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^HB-", col.name = "percent.hb")
list.samples <- lapply(X = list.samples, FUN = PercentageFeatureSet, pattern = "^RP[SL]", col.name = "percent.ribo")

head(list.samples[[2]])
rownames(list.samples[[1]])[grep("^MT-", rownames(list.samples[[1]]))]
# no mt hb or rp genes present

VlnPlot(list.samples[[3]], features = c("nCount_Spatial", "nFeature_Spatial", 
                                        "percent.mt", "percent.ribo",
                            "percent.hb"), pt.size = 0.1, ncol = 2) + NoLegend()

v1 <- VlnPlot(list.samples[[1]], features = "nCount_Spatial", pt.size = 0.1) + 
  NoLegend() + ggtitle(paste0("nCount sample ", names(list.samples)[1]))
v2 <- VlnPlot(list.samples[[2]], features = "nCount_Spatial", pt.size = 0.1) +
NoLegend() + ggtitle(paste0("nCount sample ", names(list.samples)[2]))
v3 <- VlnPlot(list.samples[[3]], features = "nCount_Spatial", pt.size = 0.1) +
NoLegend() + ggtitle(paste0("nCount sample ", names(list.samples)[3]))
v4 <- VlnPlot(list.samples[[4]], features = "nCount_Spatial", pt.size = 0.1) +
NoLegend() + ggtitle(paste0("nCount sample ", names(list.samples)[4]))
f1 <- SpatialFeaturePlot(list.samples[[1]], features = "nCount_Spatial") + 
  theme(legend.position = "right") + ggtitle(paste0(names(list.samples)[1]))
f2 <- SpatialFeaturePlot(list.samples[[2]], features = "nCount_Spatial") + 
  theme(legend.position = "right") + ggtitle(paste0(names(list.samples)[2]))
f3 <- SpatialFeaturePlot(list.samples[[3]], features = "nCount_Spatial") + 
  theme(legend.position = "right") + ggtitle(paste0(names(list.samples)[3]))
f4 <- SpatialFeaturePlot(list.samples[[4]], features = "nCount_Spatial") + 
  theme(legend.position = "right") + ggtitle(paste0(names(list.samples)[4]))

wrap_plots(v1, v2, v3, v4,
           f1, f2, f3, f4,
           ncol = 4)

SpatialFeaturePlot(list.samples[[1]], features = c("CD4", "PTPRC","CD8A", "EPCAM", "ACTA2", "FAP"))
SpatialFeaturePlot(list.samples[[2]], features = c("CD4", "PTPRC","CD8A", "EPCAM", "ACTA2", "FAP"))
SpatialFeaturePlot(list.samples[[3]], features = c("CD4", "PTPRC","CD8A", "EPCAM", "ACTA2", "FAP"))
SpatialFeaturePlot(list.samples[[4]], features = c("CD4", "PTPRC","CD8A", "EPCAM", "ACTA2", "FAP"))


## preprocess in a list ----

lapply(list.samples, dim)
list.samples <- lapply(list.samples, function(x) {
  x = x[ , x$nFeature_Spatial > 1000]
})
list.samples <- lapply(list.samples, SCTransform, assay = "Spatial", method = "poisson")
list.samples <- lapply(list.samples, RunPCA, assay = "SCT", verbose = FALSE)
list.samples <- lapply(list.samples, FindNeighbors, reduction = "pca", dims = 1:30)
list.samples <- lapply(list.samples, FindClusters, resolution = 1.2, verbose = FALSE) # increased resolution
list.samples <- lapply(list.samples, RunUMAP, reduction = "pca", dims = 1:30)

save(list.samples, file="spatial.list.normalized.RData")



## merging ----

names(list.samples)
sc.all <- merge(list.samples[[1]], list.samples[2:4], add.cell.ids = names(list.samples)) # 17943 features across 4526 samples
names(sc.all)
head(sc.all[[]])
table(sc.all$sample)
# 01_034_C1d1 01_034_C3d1 01_039_C1d1 01_039_C3d1 
#   539         755        1091        2141
sc.all$Sample_ID <- sc.all$sample
sc.all$Sample_ID <- gsub("01_034_C1d1", names(list.samples)[1], sc.all$Sample_ID)
sc.all$Sample_ID <- gsub("01_034_C3d1", names(list.samples)[2], sc.all$Sample_ID)
sc.all$Sample_ID <- gsub("01_039_C1d1", names(list.samples)[3], sc.all$Sample_ID)
sc.all$Sample_ID <- gsub("01_039_C3d1", names(list.samples)[4], sc.all$Sample_ID)
sc.all$Sample_ID <- factor(sc.all$Sample_ID, levels = names(list.samples))
table(sc.all$Sample_ID)

VlnPlot(sc.all, features = c("nCount_Spatial", "nFeature_Spatial"), 
#                             "percent.mito", "percent.hb"), 
        group.by = "Sample_ID",
        pt.size = 0.1, ncol = 2) + NoLegend() 

SpatialFeaturePlot(sc.all, features = c("nCount_Spatial", "nFeature_Spatial"))


## filtering ----

sc.all <- sc.all[, sc.all$nFeature_Spatial > 1000]

SpatialFeaturePlot(sc.all, features = c("nCount_Spatial", "nFeature_Spatial"))


# top expressed genes

C = as.matrix(sc.all@assays$Spatial@counts)
C[1:5,1:10]
dim(C)
most_expressed <- as.data.frame(rowSums(C))
head(most_expressed)
most_expressed <- most_expressed[order(most_expressed$`rowSums(C)`, decreasing = T), , drop = F]
most_expressed <- most_expressed[1:20, ,drop=F]

par(mar=c(5,10,5,5))
boxplot(t(C[rownames(most_expressed), ]), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
        main = "Top Expressed Genes")


# Filter genes

# no genes to filter

dim(sc.all)
## 17943  4200



## normalization ----

sc.all <- SCTransform(sc.all, assay = "Spatial", verbose = TRUE, method = "poisson") %>%
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

SpatialFeaturePlot(sc.all, features = c("CD4", "CD8A","EPCAM", "ACTA2"))
SpatialFeaturePlot(sc.all, features = c("CD4", "CD8A","EPCAM", "ACTA2"),
                   pt.size.factor = 2, alpha = c(0.5, 1))

DimPlot(sc.all, reduction = "umap", group.by = c("ident", "Sample_ID"))
DimPlot(sc.all, reduction = "umap", split.by = "Sample_ID")

SpatialDimPlot(sc.all)

SpatialDimPlot(sc.all, cells.highlight = CellsByIdentities(sc.all), 
               facet.highlight = TRUE,
               ncol = 5)



save(sc.all, file = "spatial.all.RData")



# consider integration



# integration -------------------------------------------------------------

rm(list=ls())

load("spatial.all.RData")

st.list <- SplitObject(sc.all, split.by = "Sample_ID")
st.list <- lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")
st.features  <- SelectIntegrationFeatures(st.list, nfeatures = 3000)
# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, 
                              anchor.features = st.features,
                              verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
sc.int <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

rm(int.anchors, st.list)
gc()

sc.int <- RunPCA(sc.int, verbose = FALSE) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

DimPlot(sc.int, reduction = "umap", group.by = c("ident", "Sample_ID"))
sc.int$Sample_ID <- factor(sc.int$Sample_ID, 
                           levels = c("P034_Pre" , "P034_Post", "P039_Pre",  "P039_Post"))
DimPlot(sc.int, reduction = "umap", split.by = "Sample_ID")

SpatialDimPlot(sc.int)


save(sc.int, file = "spatial.integrated.RData")



## spatially variable genes ----

sc.int <- FindSpatiallyVariableFeatures(sc.int, assay = 'SCT', 
                                        features = VariableFeatures(sc.int)[1:1000], 
                                        selection.method = 'markvariogram')

top.features <- head(SpatiallyVariableFeatures(sc.int, selection.method = 'markvariogram'), 10)




# annotation -----------------------------------------------------------

## integration with scRNAseq data ----

rm(list=ls())

# using individual spatial datasets already preprocessed
load("spatial.list.normalized.RData")
st.list <- list.samples

# reference dataset
load("sc.all.celltypes.RData")
Idents(sc.merged) <- sc.merged$cell_type

# select x cells per subclass
sc.merged <- subset(sc.merged, cells = WhichCells(sc.merged, downsample = 200)) # similar results using all cells
table(sc.merged$cell_type)
DefaultAssay(sc.merged) <- "SCT"


# query in a loop

st.list.anno <- list()

for(i in 1:length(st.list)){
  
  sc.query <- st.list[[i]]
  
  anchors <- FindTransferAnchors(reference = sc.merged, 
                                 query = sc.query, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = sc.merged$cell_type, 
                                    prediction.assay = TRUE,
                                    weight.reduction = sc.query[["pca"]], dims = 1:30)
  sc.query[["predictions"]] <- predictions.assay
  
  predictions.label <- TransferData(anchorset = anchors, 
                                    refdata = sc.merged$cell_type, 
                                    prediction.assay = F,
                                    weight.reduction = sc.query[["pca"]], dims = 1:30)
  sc.query <- AddMetaData(sc.query, predictions.label)
  
  st.list.anno[[i]] <- sc.query
  
  rm(sc.query, anchors, predictions.assay, predictions.label)
  
}

names(st.list.anno) <- names(st.list)

table(st.list.anno$P034_Pre$predicted.id)
# CAFs Tumor cells 
# 438          62 
table(st.list.anno$P034_Post$predicted.id)
# CAFs Tumor cells 
# 295         449 
table(st.list.anno$P039_Pre$predicted.id)
# CAFs Endothelial Tumor cells 
# 760          17          42 
table(st.list.anno$P039_Post$predicted.id)
# CAFs 
# 2137 

# no immune cells detected in any slide
# mostly CAFs and tumor cells


save(st.list.anno, file="st.list.anno.RData")



## add manual annotations ----

sc.query <- st.list.anno$P034_Pre
sc.query <- st.list.anno$P034_Post
sc.query <- st.list.anno$P039_Pre
sc.query <- st.list.anno$P039_Post

Idents(sc.query) <- sc.query$predicted.id

SpatialDimPlot(sc.query, images = "slice1",
               label = T, label.box = T, repel = T)

#DefaultAssay(sc.query) <- "Spatial"
DefaultAssay(sc.query) <- "SCT"

SpatialFeaturePlot(sc.query, c("EPCAM","ACTA2",
                               "PTPRC","CD14"), 
                   images = "slice1")

SpatialFeaturePlot(sc.query, c("CD4","CD8A","PTPRC"), images = "slice1") # immune markers
SpatialFeaturePlot(sc.query, c("EPCAM","PGR","TFF3"), images = "slice1") # tumor markers
SpatialFeaturePlot(sc.query, c("ACTA2","COL1A1","FN1"), images = "slice1") # fibroblast markers
SpatialFeaturePlot(sc.query, c("PECAM1","CD34","PROM1"), images = "slice1") # endothelial
SpatialFeaturePlot(sc.query, c("CD163","CD68","CD14"), images = "slice1") # macrophages

Idents(sc.query) <- "Undefined"

sel.cells <- WhichCells(sc.query, expression = ACTA2 > 0 & COL1A1 > 0 & FN1 > 0 & PTPRC == 0)
Idents(sc.query, cells = sel.cells) <- "CAFs"

#sel.cells <- WhichCells(sc.query, expression = EPCAM > 0 & TFF3 > 0 & ACTA2 == 0)
sel.cells <- WhichCells(sc.query, expression = EPCAM > 0  & ACTA2 == 0)
Idents(sc.query, cells = sel.cells) <- "Tumor cells"

sel.cells <- WhichCells(sc.query, expression = PTPRC > 0 & CD4 > 0 & CD8A == 0)
Idents(sc.query, cells = sel.cells) <- "CD4 T cells"

sel.cells <- WhichCells(sc.query, expression = PTPRC > 0 & CD8A > 0 & CD4 == 0)
Idents(sc.query, cells = sel.cells) <- "CD8 T cells"

sel.cells <- WhichCells(sc.query, expression = CD163 > 0 & CD14 > 0 & CD68 > 0 & PTPRC == 0)
Idents(sc.query, cells = sel.cells) <- "Macrophages"

sel.cells <- WhichCells(sc.query, expression = CD34 > 0 & PECAM1 > 0 & PROM1 > 0)
Idents(sc.query, cells = sel.cells) <- "Endothelial"

sc.query$cell_type <- Idents(sc.query)
table(sc.query$cell_type)

SpatialDimPlot(sc.query, images = "slice1",
               label = T, label.box = T, repel = T)


sc.query.1 <- sc.query
sc.query.2 <- sc.query
sc.query.3 <- sc.query
sc.query.4 <- sc.query

new.anno.list <- list(sc.query.1,
                      sc.query.2,
                      sc.query.3,
                      sc.query.4)
names(new.anno.list) <- c("P034_Pre",
                          "P034_Post",
                          "P039_Pre",
                          "P039_Post")

head(new.anno.list$P034_Pre[[]])


save(new.anno.list, file="st.list.anno.manual.RData")



## visualizations ----


lapply(new.anno.list, SpatialDimPlot, label = T, label.box = T, repel = T)

lapply(new.anno.list, SpatialFeaturePlot, c("CD4","CD8A","PTPRC")) # immune markers
lapply(new.anno.list, SpatialFeaturePlot, c("EPCAM","PGR","TFF3")) # tumor markers
lapply(new.anno.list, SpatialFeaturePlot, c("ACTA2","COL1A1","FN1")) # fibroblast markers
lapply(new.anno.list, SpatialFeaturePlot, c("PECAM1","CD34","PROM1")) # endothelial
lapply(new.anno.list, SpatialFeaturePlot, c("CD163","CD68","CD14")) # macrophages




## SCDC deconvolution (not run) ----
# gave similar results

inst = installed.packages()

if (!("xbioc" %in% rownames(inst))) {
  remotes::install_github("renozao/xbioc", dependencies = FALSE)
}
if (!("SCDC" %in% rownames(inst))) {
  remotes::install_github("meichendong/SCDC", dependencies = FALSE)
}

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))


sc.merged@active.assay = "RNA"

markers_sc <- FindAllMarkers(sc.merged, only.pos = TRUE, logfc.threshold = 0.1,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")


# select 1st sample
sc.query <- st.list$P034_Pre

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(sc.query), ]
head(markers_sc)

# Select top genes per cluster, select top by first p-value, 
# then absolute diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 + 1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))

eset_SC <- ExpressionSet(assayData = as.matrix(sc.merged@assays$RNA@counts[m_feats,
]), phenoData = AnnotatedDataFrame(sc.merged@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(sc.query@assays$Spatial@counts[m_feats,
]), phenoData = AnnotatedDataFrame(sc.query@meta.data))

deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, 
                                     sc.eset = eset_SC, 
                                     ct.varname = "cell_type",
                                     ct.sub = as.character(unique(eset_SC$cell_type)))

head(deconvolution_crc$prop.est.mvw)
test <- as.data.frame(deconvolution_crc$prop.est.mvw)
head(test)
colnames(test)
hist(test$`Tumor cells`)
table(test$`Tumor cells` > 0)


sc.query@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(sc.query@assays$SCDC@key) == 0) {
  sc.query@assays$SCDC@key = "scdc_"
}

DefaultAssay(sc.query) <- "SCDC"
head(sc.query[[]])

SpatialDimPlot(sc.query)
SpatialFeaturePlot(sc.query, 
                   features = c("Tumor cells"), 
                   pt.size.factor = 2, ncol = 2,
                   crop = TRUE)

sc.query <- FindSpatiallyVariableFeatures(sc.query, assay = "SCDC", selection.method = "markvariogram",
                                          features = rownames(sc.query), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(sc.query), 4)
SpatialPlot(object = sc.query, features = top.clusters, ncol = 2)
SpatialPlot(object = sc.query, features = top.clusters[3], label = T, ncol = 2)
VlnPlot(sc.query, group.by = "seurat_clusters", 
        features = top.clusters, pt.size = 0,
        ncol = 2)





# end ---------------------------------------------------------------------
sessionInfo()


