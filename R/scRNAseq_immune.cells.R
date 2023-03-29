
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





# end ---------------------------------------------------------------------







