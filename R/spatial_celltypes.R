
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

load("st.list.anno.manual.RData")

lapply(new.anno.list, SpatialDimPlot, label = T, label.box = T, repel = T)

lapply(new.anno.list, SpatialFeaturePlot, c("CD4","CD8A","PTPRC")) # immune markers
lapply(new.anno.list, SpatialFeaturePlot, c("EPCAM","PGR","TFF3")) # tumor markers
lapply(new.anno.list, SpatialFeaturePlot, c("ACTA2","COL1A1","FN1")) # fibroblast markers
lapply(new.anno.list, SpatialFeaturePlot, c("PECAM1","CD34","PROM1")) # endothelial
lapply(new.anno.list, SpatialFeaturePlot, c("CD163","CD68","CD14")) # macrophages


## proportions ----

output.dir <- "results/spatial/seurat/all.cells/"

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("gray","red")

load("data/objects/spatial/st.list.anno.manual.RData")

new.anno.list

sp.merged <- merge(new.anno.list$P034_Pre, c(new.anno.list$P034_Post,
                                             new.anno.list$P039_Pre,
                                             new.anno.list$P039_Post))

table(sp.merged$cell_type, sp.merged$sample)

prop.table(table(sp.merged$cell_type, sp.merged$sample))

tab1 <- table(sp.merged$sample, sp.merged$cell_type)
rownames(tab1) <- c("P34_Pre","P34_Post",
                    "P39_Pre","P39_Post")

tab34 <- tab1[1:2, ]
prop.pvalue <- c(1:ncol(tab34))
for(i in 1:ncol(tab34)){
  prop.res <- prop.test(tab34[,i], rowSums(tab34))
  prop.pvalue[i] <- prop.res$p.value
}
tab34 <- rbind(tab34, prop.pvalue)
write.csv(tab34, file = paste0(output.dir, "cluster.proportions.P34.csv"))

tab39 <- tab1[3:4, ]
prop.pvalue <- c(1:ncol(tab39))
for(i in 1:ncol(tab39)){
  prop.res <- prop.test(tab39[,i], rowSums(tab39))
  prop.pvalue[i] <- prop.res$p.value
}
tab39 <- rbind(tab39, prop.pvalue)
write.csv(tab39, file = paste0(output.dir, "cluster.proportions.P39.csv"))

tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cell_Type", "Sample","Proportion")
head(tab1)

p1 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Cell Proportions") + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_col() +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Cell Numbers") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))


png(paste0(output.dir, "Barplot.cellclass.png"))
jpeg(paste0(output.dir, "Barplot.cellclass.jpeg"), 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()



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




# other markers ----

load("st.list.anno.manual.RData")

lapply(new.anno.list, SpatialDimPlot, label = T, label.box = T, repel = T)

lapply(new.anno.list, SpatialFeaturePlot, c("CCR5")) # checkpoint marker


# EMT ----

output.dir <- "results/spatial/seurat/tumor.cells/"


library(pals)
colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("gray","red")


# PanCancer EMT signature (https://pubmed.ncbi.nlm.nih.gov/26420858/)

PanEMT <- read.csv("data/genelists/PanCancer_EMT_signature.csv")
head(PanEMT)
PanEMT$Symbol # 77

for(i in 1:length(new.anno.list)){
  sel.sample <- new.anno.list[[i]]
  sel.PanEMT <- intersect(PanEMT$Symbol, rownames(sel.sample)) # 76
  sel.sample <- AddModuleScore(sel.sample, list(sel.PanEMT), name = "PanEMT_Score")
  new.anno.list[[i]] <- sel.sample
}


lapply(new.anno.list, SpatialFeaturePlot, c("PanEMT_Score1"))
lapply(new.anno.list, VlnPlot, c("PanEMT_Score1"))


# use integrated dataset

load("data/objects/spatial/spatial.integrated.RData")

DimPlot(sc.int, group.by = "sample")

DefaultAssay(sc.int) <- "SCT"

sel.PanEMT <- intersect(PanEMT$Symbol, rownames(sc.int)) # 76
sc.int <- AddModuleScore(sc.int, list(sel.PanEMT), name = "PanEMT_Score")

VlnPlot(sc.int, "PanEMT_Score1", 
        group.by = "sample",
        cols = rep(condition.colors, 2)) +
  ggtitle("PanCancer EMT Score") + xlab("")
ggsave(paste0(output.dir, "Violins.PanEMT.all.cells.jpeg"))
ggsave(paste0(output.dir, "Violins.PanEMT.all.cells.pdf"))

library(ggpubr)
sc.int.34 <- subset(sc.int, sample == "01_034_C1d1" |
                    sample == "01_034_C3d1" )
v1 <- VlnPlot(sc.int.34, "PanEMT_Score1", 
        group.by = "sample",
        cols = condition.colors) +
  ggtitle("PanCancer EMT Score") +
  stat_compare_means() + xlab("")
v1
ggsave(paste0(output.dir, "Violins.PanEMT.P34.all.cells.jpeg"))
ggsave(paste0(output.dir, "Violins.PanEMT.P34.all.cells.pdf"))

sc.int.39 <- subset(sc.int, sample == "01_039_C1d1" |
                      sample == "01_039_C3d1" )
v2 <- VlnPlot(sc.int.39, "PanEMT_Score1", 
        group.by = "sample",
        cols = condition.colors) +
  ggtitle("PanCancer EMT Score") +
  stat_compare_means() + xlab("")
v2
ggsave(paste0(output.dir, "Violins.PanEMT.P39.all.cells.jpeg"))
ggsave(paste0(output.dir, "Violins.PanEMT.P39.all.cells.pdf"))

library(gridExtra)
grid.arrange(v1,v2, ncol=2)


# select tumor cells

sample1 <- read.csv("data/spatial/histology.labels/01-034 C1D1 Tumor.csv", row.names = 1)
sample2 <- read.csv("data/spatial/histology.labels/01-034 C3D1 Tumor.csv", row.names = 1)
sample3 <- read.csv("data/spatial/histology.labels/01-039 C1D1 Tumor.csv", row.names = 1)
sample4 <- read.csv("data/spatial/histology.labels/01-039 C3D1 Tumor.csv", row.names = 1)

sel.cells <- c(paste0("P034_Pre_", rownames(sample1)),
               paste0("P034_Post_", rownames(sample2)),
               paste0("P039_Pre_", rownames(sample3)),
               paste0("P039_Post_", rownames(sample4)))

table(sc.int$integrated_snn_res.0.8)
Idents(sc.int) <- sc.int$seurat_clusters
sc.int <- SetIdent(sc.int, sel.cells, "Tumor cells")
sc.int$histo <- Idents(sc.int)
table(sc.int$histo)
#save(sc.int , file = "data/objects/spatial/sc.int.RData")

subset <- subset(sc.int, ident = "Tumor cells")
saveRDS(subset , file = "data/objects/spatial/sc.int.tumor.cells.rds")


## only tumor cells ----

rm(list=ls())

output.dir <- "results/spatial/seurat/tumor.cells/EMT/"

condition.colors = rep(c("gray","red"), 2)

subset <- readRDS("data/objects/spatial/sc.int.tumor.cells.rds")
table(subset$Sample_ID)
DimPlot(subset, group.by = "Sample_ID")


### MsigDb Halmark ----

# from: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

library(msigdbr)
gsets_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gsets_msigdb <- gsets_msigdb[gsets_msigdb$gs_name==
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt.genes <- gsets_msigdb$gene_symbol # 204
emt.genes <- intersect(emt.genes, rownames(subset)) # 189

# score with AddModuleScore
subset <- AddModuleScore(subset, 
                        name = "EMT_score",
                        list(emt.genes))

v1 <- VlnPlot(subset, features = "EMT_score1", 
              group.by = "Sample_ID", cols = condition.colors) +
  ggtitle("EMT Score [Seurat]") + xlab("") +
#  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()


# compare to UCell:
library(UCell)
emt.genes <- list(EMT = emt.genes)
subset <- AddModuleScore_UCell(subset, features = emt.genes)
signature.names <- paste0(names(emt.genes), "_UCell")

v2 <- VlnPlot(subset, features = signature.names, 
              group.by = "Sample_ID", cols = condition.colors) +
  ggtitle("EMT Score [UCell]") + xlab("") +
#  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()

v1 + v2
grid.arrange(v1,v2, ncol = 2)
ggsave(paste0(output.dir, "Violins.EMT.Hallmark.pdf"))
ggsave(paste0(output.dir, "Violins.EMT.Hallmark.jpeg"))

SpatialFeaturePlot(subset, "EMT_UCell")

DefaultAssay(subset) <- "SCT"
SpatialFeaturePlot(subset, "EMT_UCell", crop = F,
                   min.cutoff = 0.1, max.cutoff = 0.3,
                   images = c("slice1","slice1_P034_Post.1")) +
#  ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 0.5, 1.0)) +
  ggtitle("UCell EMT HALLMARK P034-Pre")


SpatialFeaturePlot(subset, "EMT_UCell", crop = F,  images = c("slice1")) +
  ggtitle("UCell EMT HALLMARK P034-Pre")
SpatialFeaturePlot(subset, "EMT_UCell", crop = F,  images = "slice1_P034_Post.1")+
  ggtitle("UCell EMT HALLMARK P034-Post")
SpatialFeaturePlot(subset, "EMT_UCell", crop = F, images = "slice1_P039_Pre.2")+
  ggtitle("UCell EMT HALLMARK P039-Pre")
SpatialFeaturePlot(subset, "EMT_UCell", crop = F,  images = "slice1_P039_Post.3")+
  ggtitle("UCell EMT HALLMARK P039-Post")

SpatialFeaturePlot(subset, "EMT_score1", images = "slice1")
SpatialFeaturePlot(subset, "EMT_score1", images = "slice1_P034_Post.1")
SpatialFeaturePlot(subset, "EMT_score1", images = "slice1_P039_Pre.2")
SpatialFeaturePlot(subset, "EMT_score1", images = "slice1_P039_Post.3")


### PanCancer EMT signature ----

# PanCancer EMT signature (https://pubmed.ncbi.nlm.nih.gov/26420858/)

PanEMT <- read.csv("data/genelists/PanCancer_EMT_signature.csv")
head(PanEMT)
PanEMT$Symbol # 77

PanEMT <- intersect(PanEMT$Symbol, rownames(subset)) # 76

subset <- AddModuleScore(subset, list(PanEMT), name = "PanEMT_Score")

v1 <- VlnPlot(subset, features = "PanEMT_Score1", 
              group.by = "Sample_ID", cols = condition.colors) +
  ggtitle("PanCancer EMT [Seurat]") + xlab("") +
  #  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()


# compare to UCell:
library(UCell)
pan.emt.genes <- list(PanEMT = PanEMT)
subset <- AddModuleScore_UCell(subset, features = pan.emt.genes)
signature.names <- paste0(names(pan.emt.genes), "_UCell")

v2 <- VlnPlot(subset, features = signature.names, 
              group.by = "Sample_ID", cols = condition.colors) +
  ggtitle("PanCancer EMT [UCell]") + xlab("") +
  #  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))+
  RotatedAxis()

v1 + v2
grid.arrange(v1,v2, ncol = 2)
ggsave(paste0(output.dir, "Violins.EMT.PanCancer.pdf"))
ggsave(paste0(output.dir, "Violins.EMT.PanCancer.jpeg"))

SpatialFeaturePlot(subset, "PanEMT_UCell")
SpatialFeaturePlot(subset, "PanEMT_UCell", crop = F,  images = "slice1") +
  ggtitle("UCell EMT PanCancer P034-Pre")
SpatialFeaturePlot(subset, "PanEMT_UCell", crop = F,  images = "slice1_P034_Post.1")+
  ggtitle("UCell EMT PanCancer P034-Post")
SpatialFeaturePlot(subset, "PanEMT_UCell", crop = F,  images = "slice1_P039_Pre.2")+
  ggtitle("UCell EMT PanCancer P039-Pre")
SpatialFeaturePlot(subset, "PanEMT_UCell", crop = F,  images = "slice1_P039_Post.3")+
  ggtitle("UCell EMT PanCancer P039-Post")


saveRDS(subset, file = "data/objects/spatial/sc.int.tumor.cells.rds")

write.csv(subset@meta.data, file = paste0(output.dir, "tumor.metadata.csv"))





# end ---------------------------------------------------------------------
sessionInfo()


