
# spatial RNAseq data
# NP137 manuscript
# Visium for two samples, before and after therapy

# Differential expression and pathway analyses


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

set.seed(4)

setwd()

output.dir <- "results/spatial/seurat/all.cells/"


# DEGs and pathways -----------------------------------------------------------

load("data/objects/spatial/st.list.anno.manual.RData")

table(new.anno.list$P034_Pre$cell_type)
table(new.anno.list$P034_Post$cell_type)
table(new.anno.list$P039_Pre$cell_type)
table(new.anno.list$P039_Post$cell_type)

sel.celltypes <- c("Endothelial","Macrophages","Tumor cells","CAFs")

cell.list <- vector("list", length = length(sel.celltypes))
names(cell.list) <- sel.celltypes

for(i in 1:length(new.anno.list)){
  
  meta <- new.anno.list[[i]]@meta.data
  
  for(j in 1:length(cell.list)){
    temp.cells <- rownames(meta[meta$cell_type == sel.celltypes[j], ])
    temp.cells <- paste0(names(new.anno.list)[i], "_", temp.cells)
    cell.list[[j]] <- c(cell.list[[j]], temp.cells)
    rm(temp.cells)
    
  }
  
  rm(meta)

}
  
lapply(cell.list, length)


# use integrated dataset

load("data/objects/spatial/spatial.integrated.RData")

DimPlot(sc.int, group.by = "Sample_ID")

DefaultAssay(sc.int) <- "SCT"

Idents(sc.int) <- sc.int$Sample_ID

sc.int <- SetIdent(sc.int, cell.list$Endothelial, "Endothelial")
sc.int <- SetIdent(sc.int, cell.list$Macrophages, "Macrophages")
sc.int <- SetIdent(sc.int, cell.list$"Tumor cells", "Tumor cells")
sc.int <- SetIdent(sc.int, cell.list$CAFs, "CAFs")

table(Idents(sc.int))
sc.int$celltype <- Idents(sc.int)

sc.int <- PrepSCTFindMarkers(sc.int)

save(sc.int, file = "data/objects/spatial/spatial.integrated.RData")

sc.int <- subset(sc.int, ident = c("Endothelial","Tumor cells",
                                   "Macrophages","CAFs"), invert = F)

sc.int$celltype.condition <- paste(sc.int$celltype, sc.int$Sample_ID, sep="_")
table(sc.int$celltype.condition)
Idents(sc.int) <- "Sample_ID"

DefaultAssay(sc.int) <- "Spatial"


## Upregulated pathways ----

### Patient 1 ----

sc.sel <- subset(sc.int, ident = c("P034_Pre","P034_Post"), invert = F)

Idents(sc.sel) <- "celltype.condition"

listEnrichrDbs()
listEnrichrDbs()[grep("MSigDB", listEnrichrDbs()$libraryName), ]

dbs <- c("MSigDB_Hallmark_2020")

for (level in levels(factor(sc.sel$celltype))){
  try({
    ident1 <- paste0(level,"_P034_Post")
    ident2 <- paste0(level,"_P034_Pre")
    degs <- FindMarkers(sc.sel, 
                        ident.1 = ident1, 
                        ident.2=ident2, 
                        min.pct=0.5, 
                        logfc.threshold=1)
    degs <- degs[degs$p_val_adj < 0.05, ]

        write.csv(degs, file=paste0(output.dir,level, "_P034_DEGs.csv"))

        genes <- degs[degs$avg_log2FC > 0, ]
        genes <- rownames(genes)
        enriched <- enrichr(genes, dbs)

        write.csv(enriched, file=paste0(output.dir,level, "_P034_HALLMARK.csv"))

        rm(degs, genes, enriched)
    
  })
}



### Patient 2 ----

sc.sel <- subset(sc.int, ident = c("P039_Pre","P039_Post"), invert = F)

Idents(sc.sel) <- "celltype.condition"

#listEnrichrDbs()
#listEnrichrDbs()[grep("MSigDB", listEnrichrDbs()$libraryName), ]

dbs <- c("MSigDB_Hallmark_2020")

for (level in levels(factor(sc.sel$celltype))){
  try({
    ident1 <- paste0(level,"_P039_Post")
    ident2 <- paste0(level,"_P039_Pre")
    degs <- FindMarkers(sc.sel, 
                        ident.1 = ident1, 
                        ident.2=ident2, 
                        min.pct=0.5, 
                        logfc.threshold=1)
    degs <- degs[degs$p_val_adj < 0.05, ]
    
    write.csv(degs, file=paste0(output.dir,level, "_P039_DEGs.csv"))
    
    genes <- degs[degs$avg_log2FC > 0, ]
    genes <- rownames(genes)
    enriched <- enrichr(genes, dbs)
    
    write.csv(enriched, file=paste0(output.dir,level, "_P039_HALLMARK.csv"))
    
    rm(degs, genes, enriched)
    
  })
}


## Downregulated pathways ----

### Patient 1 ----

sc.sel <- subset(sc.int, ident = c("P034_Pre","P034_Post"), invert = F)

Idents(sc.sel) <- "celltype.condition"

listEnrichrDbs()
listEnrichrDbs()[grep("MSigDB", listEnrichrDbs()$libraryName), ]

dbs <- c("MSigDB_Hallmark_2020")

for (level in levels(factor(sc.sel$celltype))){
  try({
    ident1 <- paste0(level,"_P034_Post")
    ident2 <- paste0(level,"_P034_Pre")
    degs <- FindMarkers(sc.sel, 
                        ident.1 = ident1, 
                        ident.2=ident2, 
                        min.pct=0.25, 
                        logfc.threshold=0.25)
    degs <- degs[degs$p_val_adj < 0.05, ]
    
#    write.csv(degs, file=paste0(output.dir,level, "_P034_DEGs.csv"))
    
    genes <- degs[degs$avg_log2FC < 0, ]
    genes <- rownames(genes)
    enriched <- enrichr(genes, dbs)

    write.csv(enriched, file=paste0(output.dir, "pathways.down/", level, "_P034_HALLMARK.csv"))
    
    rm(degs, genes, enriched)
    
  })
}



### Patient 2 ----

sc.sel <- subset(sc.int, ident = c("P039_Pre","P039_Post"), invert = F)

Idents(sc.sel) <- "celltype.condition"

#listEnrichrDbs()
#listEnrichrDbs()[grep("MSigDB", listEnrichrDbs()$libraryName), ]

dbs <- c("MSigDB_Hallmark_2020")

for (level in levels(factor(sc.sel$celltype))){
  try({
    ident1 <- paste0(level,"_P039_Post")
    ident2 <- paste0(level,"_P039_Pre")
    degs <- FindMarkers(sc.sel, 
                        ident.1 = ident1, 
                        ident.2=ident2, 
                        min.pct=0.25, 
                        logfc.threshold=0.25)
    degs <- degs[degs$p_val_adj < 0.05, ]
    
#    write.csv(degs, file=paste0(output.dir,level, "_P039_DEGs.csv"))
    
    genes <- degs[degs$avg_log2FC < 0, ]
    genes <- rownames(genes)
    enriched <- enrichr(genes, dbs)
    
    write.csv(enriched, file=paste0(output.dir, "pathways.down/",level, "_P039_HALLMARK.csv"))
    
    rm(degs, genes, enriched)
    
  })
}




# end ---------------------------------------------------------------------
sessionInfo()


