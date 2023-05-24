
# intro -------------------------------------------------------------------

# immune profiling, NP137 manuscript

# Macrophages, CAFs



# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")

suppressPackageStartupMessages({
  library(Seurat);  
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


# macrophages ----

output.dir <- "results/scRNAseq/seurat/macrophages/"

load("data/objects/immune.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("gray","red")

Idents(immune) <- immune$cell_type
DefaultAssay(immune) <- "SCT"
DimPlot(immune)

macro <- subset(immune, idents = "Macrophages")
table(macro$ID)
# C1D1 C3D1 
# 577  216 

macro <- SCTransform(macro, vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")

DimPlot(macro, split.by = "ID")
ggsave(paste0(output.dir, "DimPlot.jpeg"))
ggsave(paste0(output.dir, "DimPlot.pdf"))

FeaturePlot(macro, c("CD163","MRC1"), split.by = NULL)
ggsave(paste0(output.dir, "FeaturePlots.jpeg"))
ggsave(paste0(output.dir, "FeaturePlots.pdf"))

FeaturePlot(macro, c("CD163","MRC1"), split.by = "ID")
VlnPlot(macro, c("CD163","MRC1"), split.by = "ID")
ggsave(paste0(output.dir, "FeaturePlots_split.jpeg"))
ggsave(paste0(output.dir, "FeaturePlots_split.pdf"))

M2 <- WhichCells(macro, 
                 expression = CD163 > 0 &
                 MRC1 > 0 &
                   CD40 == 0 &
                   CD86 == 0
                 )

DimPlot(macro, cells.highlight = M2, split.by = "ID")

Idents(macro) <- macro$cell_type
macro <- SetIdent(macro, M2, "M2")
macro <- RenameIdents(macro, "Macrophages" = "Others")

DimPlot(macro, cols = condition.colors, split.by = "ID")

macro$M_subtype <- Idents(macro)
table(macro$M_subtype, macro$ID)
#      C1D1 C3D1
# Others  481   193
# M2      96  23

# 16 vs 10%, p = 0.04649
# 
prop.test(c(96,23), c(481+96, 193+23), 
          p = NULL, alternative = "two.sided",
          correct = TRUE)

save(macro, file="data/objects/macro.RData")

prop.table(table(macro$ID, macro$M_subtype))

tab1 <- table(macro$ID, macro$M_subtype)
write.csv(tab1, file = paste0(output.dir, "M_type.proportions.csv"))
tab1 <- as.data.frame(t(tab1))
colnames(tab1) <- c("Type", "Sample","Proportion")
head(tab1)
p1 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = condition.colors) +
  ggtitle("Type Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Type Proportions") + NoLegend()
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Type)) +
  geom_col() +
  scale_fill_manual(values = condition.colors) +
  ggtitle("Type Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Type Numbers")

png(paste0(output.dir, "Barplot.M_type.png"))
jpeg(paste0(output.dir, "Barplot.M_type.jpeg"), 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()

write.csv(macro@meta.data, file = paste0(output.dir, "macro.metadata.csv"))


## score with AddModuleScore ----

# M1/M2 classification has been revisited in light of single cell data:
# https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00094-1
# new TAM subtype genelists below:
# RTM-TAMs were excluded because no specific signature is described

IFN.TAMs <- c("CASP1", "CASP4", "CCL2","CCL3","CCL4","CCL7","CCL8",
              "CD274","CD40", "CXCL2","CXCL3","CXCL9","CXCL10","CXCL11",
              "IDO1", "IFI6", "IFIT1","IFIT2","IFIT3",
              "IFITM1","IFITM3","IRF1","IRF7", "ISG15", "LAMP3",  
              "PDCD1LG2","TNFSF10","C1QA","C1QC","CD38","IL4I1","IFI44L")

Inflam.TAMs <- c("CCL2","CCL3","CCL4","CCL5","CCL20",
                 "CCL3L1","CCL3L3","CCL4L2","CCL4L4",
                 "CXCL1","CXCL2","CXCL3","CXCL5","CXCL8",
                 "G0S2","IL1B","IL1RN","IL6","INHBA","SPP1",
                 "KLF2","KLF6","NEDD9","PMAIP1","S100A8","S100A9") 

LA.TAMs <- c("ACP5", "AOPE", "APOC1", "ATF1",
             "C1QA","C1QB","C1QC", "CCL18", "CD163", 
             "CD36", "CD63", "CHI3L1", "CTSB","CTSD","CTSL",
             "F13A1", "FABP5", "FOLR2", "GPNMB", "IRF3",
             "LGALS3", "LIPA", "LPL", "MACRO",
             "MERTK", "MMP7","MMP9","MMP12", "MRC1", 
             "NR1H3", "NRF1", "NUPR1", "PLA2G7", "RNASE1", 
             "SPARC", "SPP1", "TFDP2", "TREM2", "ZEB1")

Angio.TAMs <- c("ADAM8", "AREG", "BNIP3", "CCL2","CCL4","CCL20",
                "CD163", "CD300E", "CD44", "CD55", "CEBPB", "CLEC5A",
                "CTSB", "EREG", "FCN1", "FLT1", "FN1", "HES1", 
                "IL1B", "IL1RN", "IL8", "MAF", "MIF", "NR1H3",
                "OLR1", "PPARG", "S100A8","S100A9","S100A12",
                "SERPINB2","SLC2A1","SPIC","SPP1","THBS1","TIMP1","VCAN","VEGFA")

Reg.TAMs <- c("CCL2","CD274","CD40","CD80","CD86",
              "CHIT1","CX3CR1","HLA-A","HLA-C", 
              "HLA-DQA1","HLA-DQB1", "HLA-DRA","HLA-DRB1","HLA-DRB5",
              "ICOSLG","IL10","ITGA4","LGALS9","MACRO","MRC1","TGFB2")

Prolif.TAMs <- c("CCNA2", "CDC45", "CDK1", "H2AFC", "HIST1H4C",
                 "HMGB1", "HMGN2", "MKI67", "RRM2", "STMN1", 
                 "TOP2A", "TUBA1B", "TUBB", "TYMS")


TAM.genelist <- list(IFN.TAMs=IFN.TAMs, 
                     Inflam.TAMs=Inflam.TAMs,
                     LA.TAMs=LA.TAMs,
                     Angio.TAMs=Angio.TAMs,
                     Reg.TAMs=Reg.TAMs,
                     Prolif.TAMs=Prolif.TAMs)

for(i in 1:length(TAM.genelist)){
  
  sel.list <- TAM.genelist[[i]]
  sel.list <- intersect(sel.list, rownames(macro))
  
  macro <- AddModuleScore(macro, 
                          name = names(TAM.genelist)[i],
                          list(sel.list))
  
  v1 <- VlnPlot(macro, features = paste0(names(TAM.genelist)[i], "1"), 
                group.by = "ID", cols = condition.colors) +
    xlab("") +
    stat_compare_means() +
    theme(axis.text.x = element_text(angle = 0, size = 16))
  ggsave(paste0(output.dir, "Violins.", names(TAM.genelist)[i],".jpeg"))

}

head(macro[[]])

# 3 TAM subtypes had lower scores after therapy:
# Angio.TAMs, IFN.TAMs, and Inflam.TAMs
# while 3 did not change: LA.TAMs, Prolif.TAMs, and Reg.TAMs.



# CAFs --------------------------------------------------------------------

rm(list=ls())

output.dir <- "results/scRNAseq/seurat/CAFs/"

load("data/objects/CAFs.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("gray","red")

CAFs <- SCTransform(CAFs, vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")

pdf(paste0(output.dir, "Dimplot1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot1.jpeg"), width = 500, height = 400)
DimPlot(CAFs, cols = cluster.colors)# + NoLegend()
dev.off()

pdf(paste0(output.dir, "Dimplot2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot2.jpeg"), width = 800, height = 400)
DimPlot(CAFs, cols = cluster.colors, split.by = "ID")# + NoLegend()
dev.off()

p1 <- DimPlot(CAFs, label = T, label.box = F, label.size = 5, cols = cluster.colors)
p2 <- DimPlot(CAFs, group.by = "ID", cols = cluster.colors)

pdf(paste0(output.dir, "Dimplot3.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot3.jpeg"), width = 1000, height = 400)
p1+p2
dev.off()


## score with AddModuleScore ----

# CAFs HH (don't run)

known.cafs <- read.csv("data/genelists/CAF_markers.csv")
head(known.cafs)

CAF.genelist <- list(myCAFs= toupper(known.cafs$myCAFs), 
                     iCAFs=toupper(known.cafs$iCAFs),
                     apCAFs=toupper(known.cafs$apCAFs))

##

DimPlot(CAFs, split.by = "ID")
table(CAFs$CAF.class, CAFs$ID)
FeaturePlot(CAFs, "ACTA2")
FeaturePlot(CAFs, "COL12A1")
FeaturePlot(CAFs, "CD74")
FeaturePlot(CAFs, "FBLN1", split.by = "ID")
VlnPlot(CAFs, "FBLN1", split.by = "ID")

mCAFs= c("MMP2", "DCN", "COL12A1", "FBLN1")
vCAFs= c("MCAM", "COL4A1", "COL18A1")
iCAFs1= c("MUSTN1", "TAGLN", "S100A4", "CXCL2")
iCAFs2= c("S100A8", "CXCL8", "SPLI")

CAF.genelist <- list(mCAFs= mCAFs,
                     vCAFs=vCAFs,
                     iCAFs1=iCAFs1,
                     iCAFs2=iCAFs2)

for(i in 1:length(CAF.genelist)){
  
  sel.list <- CAF.genelist[[i]]
  sel.list <- intersect(sel.list, rownames(CAFs))
  
  CAFs <- AddModuleScore(CAFs, 
                          name = names(CAF.genelist)[i],
                          list(sel.list))
  
  v1 <- VlnPlot(CAFs, features = paste0(names(CAF.genelist)[i], "1"), 
                group.by = "ID", cols = condition.colors) +
    xlab("") +
    stat_compare_means() +
    theme(axis.text.x = element_text(angle = 0, size = 16))
  ggsave(paste0(output.dir, "Violins.", names(CAF.genelist)[i],".jpeg"))
  
}

head(CAFs[[]])

#d1 <- plot_density(CAFs, "myCAFs1", reduction = "umap")
#d2 <- plot_density(CAFs, "iCAFs1", reduction = "umap")
#d3 <- plot_density(CAFs, "apCAFs1", reduction = "umap")
#d1+d2+d3

d1 <- plot_density(CAFs, "mCAFs1", reduction = "umap")
d2 <- plot_density(CAFs, "vCAFs1", reduction = "umap")
d3 <- plot_density(CAFs, "iCAFs11", reduction = "umap")
d4 <- plot_density(CAFs, "iCAFs21", reduction = "umap")

d1+d2+d3+d4

ggsave(paste0(output.dir, "DensityPlots_CAF.subtypes.jpeg"))
ggsave(paste0(output.dir, "DensityPlots_CAF.subtypes.pdf"))


# significantly less apCAFs after therapy, based on HH's CAFs

CAFs <- RenameIdents(CAFs,
                     '0'= "vCAFs",
                     '1'= "iCAFs.2",
                     '2'= "iCAFs.1",
                     '3'= "mCAFs")

pdf(paste0(output.dir, "Dimplot1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot1.jpeg"), width = 500, height = 400)
DimPlot(CAFs, cols = cluster.colors)
dev.off()

pdf(paste0(output.dir, "Dimplot2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot2.jpeg"), width = 800, height = 400)
DimPlot(CAFs, split.by = "ID", cols = cluster.colors)
dev.off()

CAFs$CAF.class <- Idents(CAFs)


write.csv(CAFs@meta.data, file = paste0(output.dir, "CAFs.metadata.csv"))

save(CAFs, file="data/objects/CAFs.RData")


## proportions ----

prop.table(table(CAFs$ID, CAFs$CAF.class))

tab1 <- table(CAFs$ID, CAFs$CAF.class)
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
  ggtitle("CAF Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("CAF Proportions") + NoLegend()
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("CAF Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("CAF Numbers")

png(paste0(output.dir, "Barplot.cellclass.png"))
jpeg(paste0(output.dir, "Barplot.cellclass.jpeg"), 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()


## cycling CAFs? ----

FeaturePlot(CAFs, "CDK1")
FeaturePlot(CAFs, "AURKA", split.by = "ID")
FeaturePlot(CAFs, "BIRC5", split.by = "ID")

caf.pre <- subset(CAFs, subset = ID == "C1D1")
caf.post <- subset(CAFs, subset = ID == "C3D1")

plot_density(CAFs, "CDK1", reduction = "umap")
plot_density(CAFs, "BUB1", reduction = "umap")

p1 <- plot_density(caf.pre, "BIRC5", reduction = "umap") +
  ggtitle("BIRC5 in C1D1")
p2 <- plot_density(caf.post, "BIRC5", reduction = "umap")+
  ggtitle("BIRC5 in C3D1")
p3 <- plot_density(caf.pre, "CDK1", reduction = "umap")+
  ggtitle("CDK1 in C1D1")
p4 <- plot_density(caf.post, "CDK1", reduction = "umap")+
  ggtitle("CDK1 in C3D1")

p1+p2+p3+p4




## ssGSEA ----

gene.sets <- getGeneSets(species = "Homo sapiens", library = "H", subcategory = NULL)
names(gene.sets)

ES <- enrichIt(obj = CAFs, 
               method = "ssGSEA", # "Ucell",
               gene.sets = gene.sets, 
               groups = 1000, cores = 4, 
               min.size = 20)
head(ES)
CAFs <- AddMetaData(CAFs, ES)
head(CAFs[[]])

ES2 <- data.frame(CAFs[[]], Idents(CAFs))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

output <- getSignificance(ES2, 
                          group = "ID", 
                          gene.sets = names(gene.sets),
                          fit = "T.test")
head(output)
output <- output[order(output$FDR, decreasing = F), ]

write.csv(output, file = paste0(output.dir, "CAFs.Hallmark.pathways.C3D1.vs.C1D1.csv"))


ridgeEnrichment(ES2, 
                gene.set = "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                group = "cluster", add.rug = TRUE)

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
ggsave(paste0(output.dir, "SplitEnrichment.IFNg.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.IFNg.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_PROTEIN_SECRETION")
ggsave(paste0(output.dir, "SplitEnrichment.secretion.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.secretion.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_KRAS_SIGNALING_DN")
ggsave(paste0(output.dir, "SplitEnrichment.KRAS.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.KRAS.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_HEDGEHOG_SIGNALING")
ggsave(paste0(output.dir, "SplitEnrichment.HEDGEHOG.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.HEDGEHOG.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_MITOTIC_SPINDLE")
ggsave(paste0(output.dir, "SplitEnrichment.MITOTIC_SPINDLE.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.MITOTIC_SPINDLE.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_G2M_CHECKPOINT")
ggsave(paste0(output.dir, "SplitEnrichment.MITOTIC_G2M.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.MITOTIC_G2M.pdf"))




#PCA <- performPCA(enriched = ES2, gene.sets = names(gene.sets), groups = c("ID", "cluster"))
#pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)

masterPCAPlot(ES2, gene.sets = names(gene.sets), 
              PCx = "PC1", PCy = "PC2", 
              top.contribution = 10)

sel.pathways <- output[order(output$FDR, decreasing = F), ]
head(sel.pathways, 20)
sel.pathways <- rownames(sel.pathways)[1:20]

dittoHeatmap(CAFs, genes = NULL, 
             metas = sel.pathways, 
             annot.by = "ID", 
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50))

masterPCAPlot(ES, gene.sets = colnames(ES)[1:ncol(ES)-1],
              PCx = "PC1", PCy = "PC2", top.contribution = 5) +
  theme(text=element_text(size=14),
        plot.title = element_text(size = 16)) +
  ggtitle("PC plot of Hallmark gene sets")

ggsave(paste0(output.dir, "PCplot.Hallmark.jpeg"))
ggsave(paste0(output.dir, "PCplot.Hallmark.pdf"))


write.csv(CAFs@meta.data, file = paste0(output.dir, "CAFs.metadata.csv"))

save(endothelial, file="data/objects/CAFs.RData")




# endothelial ----

# checking for enrichment in C3D1 vs C1D1

output.dir <- "results/scRNAseq/seurat/endothelial.cells/"

load("data/objects/endothelial.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("gray","red")

Idents(endothelial) <- endothelial$cell_type
DefaultAssay(endothelial) <- "SCT"
DimPlot(endothelial)

endothelial <- SCTransform(endothelial, vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")

DimPlot(endothelial, split.by = "ID")
ggsave(paste0(output.dir, "DimPlot.jpeg"))
ggsave(paste0(output.dir, "DimPlot.pdf"))

table(endothelial$ID)
# C1D1 C3D1 
# 382  180

write.csv(endothelial@meta.data, file = paste0(output.dir, "endothelial.metadata.csv"))

save(endothelial, file="data/objects/endothelial.RData")


## ssGSEA ----

gene.sets <- getGeneSets(species = "Homo sapiens", library = "H", subcategory = NULL)
names(gene.sets)

ES <- enrichIt(obj = endothelial, 
               method = "ssGSEA", # "Ucell",
               gene.sets = gene.sets, 
               groups = 1000, cores = 4, 
               min.size = 20)
head(ES)
endothelial <- AddMetaData(endothelial, ES)
head(endothelial[[]])

ES2 <- data.frame(endothelial[[]], Idents(endothelial))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

output <- getSignificance(ES2, 
                          group = "ID", 
                          gene.sets = names(gene.sets),
                          fit = "T.test")
head(output)
output <- output[order(output$FDR, decreasing = F), ]

write.csv(output, file = paste0(output.dir, "Endothelial.Hallmark.pathways.C3D1.vs.C1D1.csv"))


ridgeEnrichment(ES2, 
                gene.set = "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                group = "cluster", add.rug = TRUE)

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
ggsave(paste0(output.dir, "SplitEnrichment.IFNg.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.IFNg.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
ggsave(paste0(output.dir, "SplitEnrichment.IFNa.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.IFNa.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_ANGIOGENESIS")
ggsave(paste0(output.dir, "SplitEnrichment.Angiogenesis.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.Angiogenesis.pdf"))



#PCA <- performPCA(enriched = ES2, gene.sets = names(gene.sets), groups = c("ID", "cluster"))
#pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)

masterPCAPlot(ES2, gene.sets = names(gene.sets), 
              PCx = "PC1", PCy = "PC2", 
              top.contribution = 10)

sel.pathways <- output[order(output$FDR, decreasing = F), ]
head(sel.pathways, 20)
sel.pathways <- rownames(sel.pathways)[1:20]

dittoHeatmap(endothelial, genes = NULL, 
             metas = sel.pathways, 
             annot.by = "ID", 
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50))

masterPCAPlot(ES, gene.sets = colnames(ES)[1:ncol(ES)-1],
              PCx = "PC1", PCy = "PC2", top.contribution = 5) +
  theme(text=element_text(size=14),
        plot.title = element_text(size = 16)) +
  ggtitle("PC plot of Hallmark gene sets")

ggsave(paste0(output.dir, "PCplot.Hallmark.jpeg"))
ggsave(paste0(output.dir, "PCplot.Hallmark.pdf"))


write.csv(endothelial@meta.data, file = paste0(output.dir, "endothelial.metadata.csv"))

save(endothelial, file="data/objects/endothelial.RData")





# end ---------------------------------------------------------------------
sessionInfo()

