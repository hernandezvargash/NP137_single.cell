
# intro -------------------------------------------------------------------

# scRNAseq data NP137 manuscript
# differential expression in response to therapy
# one sample pre and post therapy



# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell")

suppressPackageStartupMessages({
  library(Seurat);  
  library(dplyr); library(xlsx);  library(stringr)
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals); library(scales)

})

set.seed(170323)


# cluster markers ----

load(file="data/objects/sc.all.celltypes.RData")

DimPlot(sc.merged)
DefaultAssay(sc.merged) <- "SCT"

markers <- FindAllMarkers(sc.merged,
                          only.pos = T,
                          min.pct = 0.5,
                          logfc.threshold = 1,
                          test.use = "wilcox",
                          max.cells.per.ident = 200)

head(markers)
markers <- markers[markers$p_val_adj < 0.05, ]
table(markers$cluster)

markers %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top.markers

colors3 = c("brown","orange")

DotPlot(sc.merged, features = unique(top.markers$gene), 
#        cols = colors3, 
        dot.scale = 4, split.by = NULL) +
  RotatedAxis()
ggsave("results/scRNAseq/seurat/all.cells/Dotplot.markers.jpeg",
       width = 40, height = 12, units = "cm")
ggsave("results/scRNAseq/seurat/all.cells/Dotplot.markers.pdf",
       width = 40, height = 12, units = "cm")

DotPlot(sc.merged, features = unique(top.markers$gene), 
        cols = colors3, 
        dot.scale = 4, split.by = "ID") +
  RotatedAxis()
ggsave("results/scRNAseq/seurat/all.cells/Dotplot2.markers.jpeg",
       width = 40, height = 20, units = "cm")
ggsave("results/scRNAseq/seurat/all.cells/Dotplot2.markers.pdf",
       width = 40, height = 20, units = "cm")


# differential expression ----

rm(list=ls())

load(file="data/objects/sc.all.celltypes.RData")

table(sc.merged$ID, sc.merged$cell_type)

output <- "results/scRNAseq/seurat/all.cells/DEGs."
sc.merged$celltype.condition <- paste(sc.merged$cell_type, sc.merged$ID, sep="_")
Idents(sc.merged) <- "celltype.condition"

for (level in levels(factor(sc.merged$cell_type))){
  try({
    ident1 <- paste0(level,"_C3D1")
    ident2 <- paste0(level,"_C1D1")
    degs <- FindMarkers(sc.merged, 
                        ident.1 = ident1, 
                        ident.2=ident2, 
                        min.pct=0.5, 
                        logfc.threshold=0.5)
    write.csv(degs, file=paste0(output,level,".csv"))
  })
}



# end ---------------------------------------------------------------------
sessionInfo()






