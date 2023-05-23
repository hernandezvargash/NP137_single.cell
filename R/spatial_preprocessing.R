
# spatial scRNAseq data
# NP137 manuscript
# Visium for two samples, before and after therapy

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




# integrate patient P034 -------------------------------------------------------------

rm(list=ls())

list.files()

samples <- list.dirs(path = "data/spatial/", full.names = F, recursive = F)

data.dir <- list.dirs(path = "data/spatial", full.names = T, recursive = F)

list.samples <- list()
for(i in 1:2){
  list.samples[[i]] <- Load10X_Spatial(data.dir = paste0(data.dir[i], "/outs/"))
  list.samples[[i]]$sample <- samples[[i]]
}

names(list.samples) <- c("P034_Pre","P034_Post")
lapply(list.samples, dim)

list.samples <- lapply(list.samples, SCTransform, assay = "Spatial", method = "poisson")
st.features  <- SelectIntegrationFeatures(list.samples, nfeatures = 3000)
# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB
st.features = SelectIntegrationFeatures(list.samples, nfeatures = 3000, verbose = FALSE)
list.samples <- PrepSCTIntegration(object.list = list.samples, 
                              anchor.features = st.features,
                              verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = list.samples, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
sc.int <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                        verbose = FALSE)

sc.int <- RunPCA(sc.int, verbose = FALSE) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

head(sc.int[[]])

Idents(sc.int) <- sc.int$sample

DimPlot(sc.int, reduction = "umap")#, group.by = c("ident", "sample"))

SpatialDimPlot(sc.int)
SpatialFeaturePlot(sc.int, "EPCAM")

saveRDS(sc.int , file = "data/objects/spatial/integrated_P034.rds")



# subset in tumor cells

sample1 <- read.csv("data/spatial/histology.labels/01-034 C1D1 Tumor.csv", row.names = 1)
sample2 <- read.csv("data/spatial/histology.labels/01-034 C3D1 Tumor.csv", row.names = 1)

sel.cells <- c(paste0(rownames(sample1),"_1"),
               paste0(rownames(sample2),"_2"))

table(sc.int$integrated_snn_res.0.8)
Idents(sc.int) <- sc.int$seurat_clusters
sc.int <- SetIdent(sc.int, sel.cells, "Tumor cells")
sc.int$histo <- Idents(sc.int)
table(sc.int$histo, sc.int$sample)
#saveRDS(sc.int , file = "data/objects/spatial/integrated_P034.rds")

subset <- subset(sc.int, ident = "Tumor cells")
saveRDS(subset , file = "data/objects/spatial/integrated_P034_tumor.cells.rds")

### MsigDb Halmark ----

# from: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

library(msigdbr)
gsets_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gsets_msigdb <- gsets_msigdb[gsets_msigdb$gs_name==
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt.genes <- gsets_msigdb$gene_symbol # 204
emt.genes <- intersect(emt.genes, rownames(subset)) # 144

# score with AddModuleScore
DefaultAssay(subset) <- "SCT"
subset <- AddModuleScore(subset, 
                         name = "EMT_score",
                         list(emt.genes))

condition.colors = rep(c("gray","red"), 2)

library(ggpubr)
v1 <- VlnPlot(subset, features = "EMT_score1", 
              group.by = "sample", cols = condition.colors) +
  ggtitle("EMT Score [Seurat]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()


# compare to UCell:
library(UCell)
emt.genes <- list(EMT = emt.genes)
subset <- AddModuleScore_UCell(subset, features = emt.genes)
signature.names <- paste0(names(emt.genes), "_UCell")

v2 <- VlnPlot(subset, features = "EMT_UCell", 
              group.by = "sample", cols = condition.colors) +
  ggtitle("EMT Score [UCell]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()

v1 + v2
library(gridExtra)
grid.arrange(v1,v2, ncol = 2)
ggsave(paste0(output.dir, "Violins.EMT.Hallmark.pdf"))
ggsave(paste0(output.dir, "Violins.EMT.Hallmark.jpeg"))


saveRDS(subset , file = "data/objects/spatial/integrated_P034_tumor.cells.rds")
subset <- readRDS("data/objects/spatial/integrated_P034_tumor.cells.rds")


SpatialFeaturePlot(subset, "EMT_UCell")


DefaultAssay(subset) <- "Spatial"

p <- SpatialFeaturePlot(subset, "EMT_UCell", crop = F)

p + ggplot2::scale_fill_continuous(limits = c(0.0,10.0),, breaks = c(0.0, 5.0, 1.0))



DefaultAssay(subset) <- "integrated"

FeaturePlot(subset, "EMT_UCell")
FeaturePlot(subset, "EMT_UCell", split.by = "sample")
VlnPlot(subset, "EMT_UCell", split.by = "sample")

subset$sample <- factor(subset$sample)

SpatialFeaturePlotSingle<- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj$EMT_UCell)
  maximal<- max(obj$EMT_UCell)
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    obj2 <- obj[, subset_cells]
    p<- SpatialFeaturePlot(obj2, features = feature, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }

  return(ps)
}

p_list<- SpatialFeaturePlotSingle(subset, feature= "EMT_UCell", 
                           metadata_column = "sample")

layout1<-"
AB
#DE
"
wrap_plots(p_list ,guides = 'collect', design = layout1)

p_list$'01_034_C1d1' + p_list$'01_034_C3d1'


minimal<- min(subset$EMT_UCell)
maximal<- max(subset$EMT_UCell)

SpatialFeaturePlot(subset, "EMT_UCell", crop = F) +
  scale_color_viridis_c(limits=c(minimal, maximal), direction = 1)

library('STutility')

FeatureOverlay(subset, features = c("LOXL1"), sampleids = c(1,2), 
               cols = c("dodgerblue2", "seagreen3", "khaki2", "darkgoldenrod2", "brown4"), 
               ncol = 2, pt.size = 4, value.scale = "all")



SpatialFeaturePlot(subset, "EMT_UCell", crop = F,  images = c("slice1")) +
  ggtitle("UCell EMT HALLMARK P034-Pre")
SpatialFeaturePlot(subset, "EMT_UCell", crop = F,  images = "slice1_P034_Post.1")+
  ggtitle("UCell EMT HALLMARK P034-Post")


save(sc.int, file = "spatial.integrated.RData")



## spatially variable genes ----

sc.int <- FindSpatiallyVariableFeatures(sc.int, assay = 'SCT', 
                                        features = VariableFeatures(sc.int)[1:1000], 
                                        selection.method = 'markvariogram')

top.features <- head(SpatiallyVariableFeatures(sc.int, selection.method = 'markvariogram'), 10)




# end ---------------------------------------------------------------------
sessionInfo()


