
# intro -------------------------------------------------------------------

# immune profiling with scRNAseq data


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")

suppressPackageStartupMessages({
  library(Seurat);  library(SeuratDisk)
  library(dplyr); library(xlsx);  library(stringr)
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals); library(scales)
  library(msigdbr); library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(biomaRt)
  library(Scillus); library(Nebulosa)
  library(escape); library(dittoSeq)
  library(UCell); library(AUCell)
  
})

set.seed(170323)

output.dir <- "results/scRNAseq/seurat/immune.cells/"


# preprocessing ---------------------------------------------

load("data/objects/immune.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:16]
condition.colors = c("brown","orange")

immune <- SCTransform(immune, vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.1) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")


DimPlot(immune, group.by = "cell_type", split.by = "ID")

DimPlot(immune, group.by = "Scibet_1")#, pt.size = 1)
DimPlot(immune, group.by = "Scibet_2")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_main")#, pt.size = 1)
DimPlot(immune, group.by = "SingleR_fine", label = T, repel = T) + NoLegend() #, pt.size = 1)
DimPlot(immune, group.by = "ScType")#, pt.size = 1)
DimPlot(immune, group.by = "cell_type")#, pt.size = 1)



# T cells -------------------------------------------------

head(immune[[]])

Idents(immune) <- immune$cell_type
table(immune$cell_type)

Tcells <- subset(immune, idents = c('CD8+ T cells',
                                    'CD4+ T cells',
                                    'NK cells'), invert = F)
table(Tcells$ID, Tcells$cell_type)

# increased resolution
Tcells <- SCTransform(Tcells, vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")

DimPlot(Tcells)

p1 <- DimPlot(Tcells, label = F, 
              group.by = "cell_type",
              label.box = F, label.size = 5)
p2 <- DimPlot(Tcells, group.by = "ID", cols = condition.colors)
p1+p2

DimPlot(Tcells, 
        group.by = "cell_type",
        split.by = "ID")#, pt.size = 1)

FeaturePlot(Tcells, c("KLRB1","CD8A","CD3D","PTPRC"))
rownames(Tcells[grep("TRG", rownames(Tcells))])
FeaturePlot(Tcells, c("TRGC1","TRGC2","TRGV2","TRGV3"))
FeaturePlot(Tcells, c("FOXP3","IL2RA","CTLA4","CD4"))


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
                                   split.by = "ID")
head(querydata[[]])
table(querydata$functional.cluster)

p1 <- DimPlot(querydata, label = T, label.box = F, label.size = 3)
p2 <- DimPlot(querydata, group.by = "cell_type", label = T) + NoLegend()
p3 <- DimPlot(querydata, group.by = "functional.cluster", label = T, repel = T) + NoLegend()

p2+p3

DimPlot(querydata, group.by = "functional.cluster", label = F, split.by = "sample")

head(querydata[[]])

# ProjecTILs wrongly re-labels NK cells,
# likely because this cell type is not present in the reference
# However, it identifies Tregs and Tex cells


## CellTypist ----

# save for CellTypist:

Tcells
DefaultAssay(Tcells) <- "RNA"
Tcells.log <- NormalizeData(Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
SaveH5Seurat(Tcells.log, filename = "data/objects/Tcells.h5Seurat")
Convert("data/objects/Tcells.h5Seurat", dest = "h5ad")


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



## new categories ----

#load("Tcells.RData")

Tcells <- querydata
rm(querydata)
head(Tcells[[]])

# automatic classifications

DimPlot(Tcells, label = T)
DimPlot(Tcells, group.by = "cell_type") # sctype
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


# only Tregs are consistently identified by several methods
# and confirmed manually as cluster 3

Tregs <- WhichCells(Tcells, idents = '3')

Idents(Tcells) <- Tcells$cell_type
Tcells <- SetIdent(Tcells, cells = Tregs, "Tregs")

pdf(paste0(output.dir, "Dimplot.Tcells.1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot.Tcells.1.jpeg"), width = 500, height = 400)
DimPlot(Tcells)
dev.off()

pdf(paste0(output.dir, "Dimplot.Tcells.2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot.Tcells.2.jpeg"), width = 800, height = 400)
DimPlot(Tcells, split.by = "ID")
dev.off()


Tcells$final.cats <- Idents(Tcells)

save(Tcells, file="data/objects/Tcells.RData")

write.csv(Tcells@meta.data, file = paste0(output.dir, "Tcells.metadata.csv"))




# final immune cell subtypes ---------------------------------------------------

#load("data/objects/immune.RData")

Idents(immune) <- immune$cell_type

immune <- SetIdent(immune, cells = Tregs, "Tregs")


pdf(paste0(output.dir, "Dimplot1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot1.jpeg"), width = 500, height = 400)
DimPlot(immune, cols = cluster.colors)# + NoLegend()
dev.off()

pdf(paste0(output.dir, "Dimplot2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot2.jpeg"), width = 800, height = 400)
DimPlot(immune, cols = cluster.colors, split.by = "ID")# + NoLegend()
dev.off()

p1 <- DimPlot(immune, label = T, label.box = F, label.size = 5, cols = cluster.colors)
p2 <- DimPlot(immune, group.by = "ID", cols = condition.colors)

pdf(paste0(output.dir, "Dimplot3.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot3.jpeg"), width = 1000, height = 400)
p1+p2
dev.off()

immune$final.cats <- Idents(immune)



## proportions ----

prop.table(table(immune$ID, immune$final.cats))

tab1 <- table(immune$ID, immune$final.cats)
prop.pvalue <- c(1:ncol(tab1))
for(i in 1:ncol(tab1)){
  prop.res <- prop.test(tab1[,i], rowSums(tab1))
  prop.pvalue[i] <- prop.res$p.value
}
tab2 <- rbind(tab1, prop.pvalue)
write.csv(tab2, file = paste0(output.dir, "cluster.proportions.csv"))

tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Cluster", "Sample","Proportion")
head(tab1)
p1 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell Type Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Cell Type Proportions") + NoLegend()
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell Type Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Cell Type Numbers")

png(paste0(output.dir, "Barplot.cellclass.png"))
jpeg(paste0(output.dir, "Barplot.cellclass.jpeg"), 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()



write.csv(immune@meta.data, file = paste0(output.dir, "immune.metadata.csv"))


save(immune, file="data/objects/immune.RData")




# checkpoint genes --------------------------------------------------------

rm(list=ls())

output.dir <- "results/scRNAseq/seurat/all.cells/"

load(file="data/objects/sc.all.celltypes.RData")
head(sc.merged[[]])
DimPlot(sc.merged, label = T)

known.markers <- c("CXCL9","CCR5","CXCL13",
                   "PDCD1","CD274",
                   "CTLA4","CD80","CD86")

jpeg(paste0(output.dir,"checkpoint.jpeg"), width = 1200, height = 900, quality = 100)
FeaturePlot(sc.merged, known.markers, ncol = 3)
dev.off()

jpeg(paste0(output.dir,"checkpoint_density.jpeg"), width = 1200, height = 900, quality = 100)
plot_density(sc.merged, reduction = "umap", known.markers)
dev.off()


i = 2
VlnPlot(sc.merged, known.markers[i], 
        group.by = "ID", ncol = 1) +
  stat_compare_means() + xlab("")
ggsave(paste0(output.dir, "Violin.CCR5_all.cells.jpeg"))
ggsave(paste0(output.dir, "Violin.CCR5_all.cells.pdf"))

FeaturePlot(sc.merged, known.markers[i], 
            order = T, pt.size = 0.6,
            split.by = "ID")
ggsave(paste0(output.dir, "CCR5_all.cells.jpeg"))
ggsave(paste0(output.dir, "CCR5_all.cells.pdf"))


VlnPlot(sc.merged, known.markers[i],
        group.by = "cell_type",
        split.by = "ID", ncol = 1)


cd8 <- subset(sc.merged, ident = "CD8+ T cells")
table(cd8$ID)

#DefaultAssay(cd8) <- "SCT"
VlnPlot(cd8, known.markers[i], 
        group.by = "ID", ncol = 1) +
  stat_compare_means() + xlab("")
ggsave(paste0(output.dir, "Violin.CCR5_CD8Tcells.jpeg"))
ggsave(paste0(output.dir, "Violin.CCR5_CD8Tcells.pdf"))




# end ---------------------------------------------------------------------







