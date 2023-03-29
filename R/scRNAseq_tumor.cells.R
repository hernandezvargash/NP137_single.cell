
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







