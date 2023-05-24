
# intro -------------------------------------------------------------------

# subset in tumor cells


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd()

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

output.dir <- "results/scRNAseq/seurat/tumor.cells/"

  
# epithelial preprocessing ----

load("data/objects/tumor.RData")

colors <- as.vector(alphabet(n=24))
cluster.colors = colors[1:12]
condition.colors = c("brown","orange")

tumor <- SCTransform(tumor, vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(dims = 1:30, umap.method = "umap-learn", metric = "correlation")

pdf(paste0(output.dir, "Dimplot1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot1.jpeg"), width = 500, height = 400)
DimPlot(tumor, cols = cluster.colors)# + NoLegend()
dev.off()

pdf(paste0(output.dir, "Dimplot2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot2.jpeg"), width = 800, height = 400)
DimPlot(tumor, cols = cluster.colors, split.by = "ID")# + NoLegend()
dev.off()

p1 <- DimPlot(tumor, label = T, label.box = F, label.size = 5, cols = cluster.colors)
p2 <- DimPlot(tumor, group.by = "ID", cols = cluster.colors)

pdf(paste0(output.dir, "Dimplot3.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot3.jpeg"), width = 1000, height = 400)
p1+p2
dev.off()


save(tumor, file="data/objects/tumor.RData")


## proportions ----

prop.table(table(tumor$ID, tumor$seurat_clusters))

tab1 <- table(tumor$ID, tumor$seurat_clusters)
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
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("") + 
  ylab("Cluster Proportions") + NoLegend()
p2 <- ggplot(tab1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Numbers") +
  theme_ipsum() +
  xlab("") +
  ylab("Cluster Numbers")

png(paste0(output.dir, "Barplot.cellclass.png"))
jpeg(paste0(output.dir, "Barplot.cellclass.jpeg"), 
     width = 500, height = 350, quality = 100)
p1 + p2
dev.off()



# ssGSEA ------------------------------------------------------------------


#load("data/objects/tumor.RData")

gene.sets <- getGeneSets(species = "Homo sapiens", library = "H", subcategory = NULL)

ES <- enrichIt(obj = tumor, 
               method = "ssGSEA", # "Ucell",
               gene.sets = gene.sets, 
               groups = 1000, cores = 4, 
               min.size = 20)
save(ES, file = "data/objects/ES_tumor.cells.RData")
head(ES)
tumor <- AddMetaData(tumor, ES)
head(tumor[[]])

ES2 <- data.frame(tumor[[]], Idents(tumor))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"


ridgeEnrichment(ES2, 
                gene.set = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                group = "cluster", add.rug = TRUE)
ggsave(paste0(output.dir, "Ridgeplot.EMT.jpeg"))
ggsave(paste0(output.dir, "Ridgeplot.EMT.pdf"))

splitEnrichment(ES2, split = "ID", gene.set = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
ggsave(paste0(output.dir, "SplitEnrichment.EMT.jpeg"))
ggsave(paste0(output.dir, "SplitEnrichment.EMT.pdf"))

PCA <- performPCA(enriched = ES2, gene.sets = names(gene.sets), groups = c("ID", "cluster"))
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)

masterPCAPlot(ES2, gene.sets = names(gene.sets), 
              PCx = "PC1", PCy = "PC2", 
              top.contribution = 10)

output <- getSignificance(ES2, 
                          group = "ID", 
                          gene.sets = names(gene.sets),
                          fit = "T.test")
head(output)
write.csv(output, file = paste0(output.dir, "Hallmark.pathways.C3D1.vs.C1D1.csv"))

sel.pathways <- output[order(output$FDR, decreasing = F), ]
head(sel.pathways, 20)
sel.pathways <- rownames(sel.pathways)[1:20]

dittoHeatmap(tumor, genes = NULL, 
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



# EMT scores ----

## MsigDb Hallmark ----

# from: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

# related signatures from MsigDb:
gsets_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gsets_msigdb <- gsets_msigdb[gsets_msigdb$gs_name==
                               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt.genes <- gsets_msigdb$gene_symbol # 204
emt.genes <- intersect(emt.genes, rownames(tumor)) # 196

# score with AddModuleScore
tumor <- AddModuleScore(tumor, 
                             name = "EMT_score",
                             list(emt.genes))

v1 <- VlnPlot(tumor, features = "EMT_score1", 
              group.by = "ID", cols = condition.colors) +
  ggtitle("EMT Score [Seurat]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))


# compare to UCell:
emt.genes <- list(EMT = emt.genes)
tumor <- AddModuleScore_UCell(tumor, features = emt.genes)
signature.names <- paste0(names(emt.genes), "_UCell")

v2 <- VlnPlot(tumor, features = "EMT_UCell", 
              group.by = "ID", cols = condition.colors) +
  ggtitle("EMT Score [UCell]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))


# compare with AUCell
cells_AUC_EMT <- AUCell_run(tumor@assays$RNA@counts, emt.genes)
tumor[["EMT_cells"]] <- cells_AUC_EMT@assays@data@listData[["AUC"]][1,]

v3 <- VlnPlot(tumor,features = "EMT_cells", 
              group.by = "ID", cols = condition.colors) +
  ggtitle("EMT Score [AUCell]") + xlab("") +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))


pdf(paste0(output.dir, "Violins.EMT.Hallmark.pdf"),
    width = 10, height = 5)
jpeg(paste0(output.dir, "Violins.EMT.Hallmark.jpeg"),
    width = 1000, height = 400)
grid.arrange(v1,v2,v3, ncol = 3)
dev.off()


VlnPlot(tumor,features = "EMT_cells", 
        split.by = "ID", cols = condition.colors,
        group.by = "seurat_clusters") +
  ggtitle("EMT Score [AUCell]") + 
  xlab("") +
  #  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 0, size = 16))
ggsave(paste0(output.dir, "Violins.by.cluster.EMT.Hallmark.pdf"))
ggsave(paste0(output.dir, "Violins.by.cluster.EMT.Hallmark.jpeg"))
              
head(tumor[[]])
d1 <- plot_density(tumor, "EMT_UCell", reduction = "umap") +
  ggtitle("EMT Score [Seurat]")
d2 <- plot_density(tumor, "EMT_score1", reduction = "umap") +
  ggtitle("EMT Score [UCell]")
d3 <- plot_density(tumor, "EMT_cells", reduction = "umap") +
  ggtitle("EMT Score [AUCell]")

d1 + d2 + d3

ggsave(paste0(output.dir, "DensityPlots.EMT.Hallmark.pdf"),
       width = 12, height = 4)
ggsave(paste0(output.dir, "DensityPlots.EMT.Hallmark.jpeg"),
       width = 12, height = 4)



## by cluster

Idents(tumor) <- tumor$seurat_clusters
DimPlot(tumor, cols = cluster.colors)# + NoLegend()
VlnPlot(tumor, "EMT_UCell", split.by = "ID", cols = condition.colors)
VlnPlot(tumor, "EMT_UCell", cols = cluster.colors)

ridgeEnrichment(ES2, 
                gene.set = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                group = "cluster", add.rug = TRUE)

library(SCpubr)

SCpubr::do_RidgePlot(sample = tumor, colors.use = cluster.colors, flip = F,
                     font.size = 12, split.by = "ID",
                     plot.title = "Hallmark EMT Score [UCell]",
                          feature = "EMT_UCell")

ES2 <- data.frame(tumor[[]], Idents(tumor))
head(ES2)
colnames(ES2)[ncol(ES2)] <- "cluster"

show_col(cluster.colors)
ridgeEnrichment(ES2, colors = cluster.colors[1:9],
                gene.set = "EMT_UCell", 
                group = "cluster", add.rug = TRUE)

splitEnrichment(ES2, split = "ID", gene.set = "EMT_UCell", colors = condition.colors)



## PanCancer EMT signature ----

# PanCancer EMT signature (https://pubmed.ncbi.nlm.nih.gov/26420858/)

PanEMT <- read.csv("data/genelists/PanCancer_EMT_signature.csv")
head(PanEMT)
PanEMT$Symbol # 77

PanEMT <- intersect(PanEMT$Symbol, rownames(tumor)) # 76

tumor <- AddModuleScore(tumor, list(PanEMT), name = "PanEMT_Score")

FeaturePlot(tumor, "PanEMT_Score1", 
            order = T, pt.size = 0.5)
plot_density(tumor, features = "PanEMT_Score1", 
             reduction = "umap",
             joint = F, pal = "plasma") + 
  ggtitle("PanCancer EMT Score")

VlnPlot(tumor, "PanEMT_Score1", cols = cluster.colors) +
  ggtitle("PanEMT Score")
ggsave(paste0(output.dir, "Violins.by.cluster.PanEMT.jpeg"))
ggsave(paste0(output.dir, "Violins.by.cluster.PanEMT.pdf"))

VlnPlot(tumor, "PanEMT_Score1", 
        split.by = "ID",
        cols = condition.colors) +
  ggtitle("PanEMT Score")
ggsave(paste0(output.dir, "Violins.by.cluster.and.condition.PanEMT.jpeg"))
ggsave(paste0(output.dir, "Violins.by.cluster.and.condition.PanEMT.pdf"))

VlnPlot(tumor, "PanEMT_Score1", group.by = "ID",
        cols = condition.colors) +
  stat_compare_means() +
  ggtitle("PanEMT Score")
ggsave(paste0(output.dir, "Violins.PanEMT.pdf"),
       width = 4, height = 5)
ggsave(paste0(output.dir, "Violins.PanEMT.jpeg"),
       width = 4, height = 5)


# clusters 5, 6, and 7 display the highest EMT scores
# but they are not present in C3D1 sample
# cluster 2 has the next highest levels
# and is present in both samples

new.idents <- c(rep("EMT.low", 2), "EMT.high",
                rep("EMT.low", 2),rep("EMT.high", 3), 
                "EMT.low")
names(new.idents) <- levels(tumor)
tumor <- RenameIdents(tumor, new.idents)

pdf(paste0(output.dir, "Dimplot.EMT1.pdf"), width = 8, height = 6)
jpeg(paste0(output.dir, "Dimplot.EMT1.jpeg"), width = 500, height = 400)
DimPlot(tumor)
dev.off()

pdf(paste0(output.dir, "Dimplot.EMT2.pdf"), width = 12, height = 6)
jpeg(paste0(output.dir, "Dimplot.EMT2.jpeg"), width = 800, height = 400)
DimPlot(tumor, split.by = "ID")
dev.off()

tumor$EMT.class <- Idents(tumor)

save(tumor, file="data/objects/tumor.RData")

write.csv(tumor@meta.data, file = paste0(output.dir, "tumor.metadata.csv"))



# end ---------------------------------------------------------------------
sessionInfo()






