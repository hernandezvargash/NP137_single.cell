
# immune profiling with scRNAseq data
# revision to NP137 manuscript
# 


# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Lab/TZones/DKO_model/integration_exp2_exp3/")

suppressPackageStartupMessages({
  library(Seurat);  library(sctransform);  library(SeuratDisk); library(SeuratWrappers); library(SingleCellExperiment)
  library(dplyr); library(xlsx);  library(stringr)
  library(ggplot2); library(cowplot); library(viridis); library(gridExtra); library(RColorBrewer); library(patchwork); library(hrbrthemes); library(pals)
  library(SingleR); library(scRNAseq); library(scibet)
  library(slingshot, quietly = T)
  library(mclust, quietly = T); library(celldex); library(KernSmooth); library(dendextend); library(circlize)
  library(enrichR); library(clusterProfiler); library(enrichplot); library(DOSE); library(pathview); library(ReactomePA)
  library(ggpubr); library(eulerr); library(alluvial)
  library(bumphunter)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(biomaRt)
  library(Scillus)
  library(SingleCellSignalR)
  library(Tempora)
  library(infercnv)
  library(SCENIC)
  library(snow)
  
})

set.seed(170323)



# subsetting after re-integration ---------------------------------------------

load(file="cc.int.sct.RData")

DefaultAssay(cc.int.sct) <- "integrated"
head(cc.int.sct[[]])
Idents(cc.int.sct) <- cc.int.sct$sample
levels(Idents(cc.int.sct))

sc.sub <- subset(cc.int.sct, idents = c("DKO_4_Lympho","DKO_8_Lympho","WT_4_Lympho","WT_8_Lympho"))
table(sc.sub$sample)
# DKO_4_Lympho DKO_8_Lympho  WT_4_Lympho  WT_8_Lympho 
# 4465          507         2488         2900 

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 100, dims = 1:30)
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "sample")

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)

Idents(sc.sub)


saveRDS(sc.sub, file = "immune.rds")




# markers -----------------------------------------------------------------

sc.sub <- readRDS("immune.rds")

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "immune.marker.reclustering.csv")
write.xlsx(all.markers, file="immune.marker.reclustering.xlsx")

all.markers <- read.csv(file = "immune.marker.reclustering.csv", row.names = 1)

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample"),
             anno_colors = list(
               alphabet2(24),
               rainbow(4)),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 5)




# automated identification ------------------------------------------------

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])

# ScType
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# get cell-type-specific gene sets from our in-built database (DB)
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare(db_, tissue)
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max = sctype_score(scRNAseqData = sc.sub[["integrated"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(sc.sub@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sc.sub@meta.data[sc.sub@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc.sub@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


# plot

sc.sub@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc.sub@meta.data$customclassif[sc.sub@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sc.sub, reduction = "umap", label = TRUE, label.size = 6, label.box = T, repel = TRUE, group.by = 'customclassif')        
DimPlot(sc.sub, reduction = "umap", label = F, repel = TRUE, group.by = 'customclassif', split.by = "sample")


# bubble plot

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(sc.sub, reduction = "umap", label = TRUE, repel = TRUE, cols = alphabet(23)), gggr, cols = 2)

saveRDS(sc.sub, file = "immune.rds")



# cell type proportions ---------------------------------------------------

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])

tab1 <- table(sc.sub$sample, sc.sub$seurat_clusters)
write.csv(tab1, file = "immune.seurat.recluster.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:24]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cell_type", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Seurat Cluster Proportions") +
  theme_ipsum() +
  xlab("")


tab2 <- table(sc.sub$sample, sc.sub$customclassif)
write.csv(tab2, file = "immune.predicted.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("New_clusters", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = New_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Predicted Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2



# filter out contaminants -------------------------------------------------

rm(list=ls())

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])
Idents(sc.sub)

DimPlot(sc.sub, label = T, label.size = 6)

FeaturePlot(sc.sub, c("Muc13","Cdh5","Krt14","Krt5","Col4a1","Col6a1","Bmp4"))
VlnPlot(sc.sub, c("Muc13","Cdh5","Krt14","Krt5","Col4a1","Col6a1","Bmp4"))

# cluster 8 : Muc13 expressing cells
# cluster 9 : Stromal cells Bmp4
# cluster 18 : Epi Cdh5
# Cluster 17 : Epi Krt14, Krt5
# Cluster 9, 20, 21, 18 : Stromal Col4a1
# Cluster 9, 20, 21 : Stromal Col6a1

sc.sub <- subset(sc.sub, idents = c('8','9','18','17','20','21'), invert = T)
table(sc.sub$sample)

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 100, dims = 1:30)
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "customclassif", label = T)

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)


saveRDS(sc.sub, file = "immune.2.rds")



## cell type proportions ---------------------------------------------------

sc.sub <- readRDS("immune.2.rds")

tab1 <- table(sc.sub$sample, Idents(sc.sub))
tab1
prop.table(tab1)
write.csv(tab1, file = "immune.ident.cluster.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cluster", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("")

tab2 <- table(sc.sub$sample, sc.sub$customclassif)
tab2
prop.table(tab2)
write.csv(tab2, file = "immune.ident.cluster.proportions.2.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("Cell_Type", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell_Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2


## markers -----------------------------------------------------------------

sc.sub <- readRDS("immune.2.rds")
DefaultAssay(sc.sub) <- "integrated"

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 1127 1124 1334 1050 1166 1446 1262 1103 1142 1335 1145 1116 1155 1408 1023  855  943  487  373   
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "immune.markers.reclustering.csv")
write.xlsx(all.markers, file="immune.markers.reclustering.xlsx")

all.markers <- read.csv(file = "immune.markers.reclustering.csv", row.names = 1)

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample","customclassif","functional.cluster"),
             anno_colors = list(
               alphabet2(20),
               rainbow(14),
               alphabet(14),
               rev(brewer.pal(n = 11, name = "Paired"))),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 4)




# T cells -------------------------------------------------

rm(list=ls())

sc.sub <- readRDS("immune.2.rds")
head(sc.sub[[]])
Idents(sc.sub)

Idents(sc.sub) <- sc.sub$customclassif
table(sc.sub$customclassif, sc.sub$sample)

DimPlot(sc.sub, label = T, label.size = 6)

sc.sub <- subset(sc.sub, idents = c('Effector CD4+ T cells'), invert = F)
table(sc.sub$sample)
# DKO_4_Lympho DKO_8_Lympho  WT_4_Lympho  WT_8_Lympho 
# 990           65          264          246 

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 20, dims = 1:30) # reduced k.weight because of low number of cells in one sample
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = T, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "customclassif", label = F)

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)


saveRDS(sc.sub, file = "Tcells.rds")



## ProjecTILs ----

# https://carmonalab.github.io/ProjecTILs_CaseStudies/
# https://github.com/carmonalab/ProjecTILs
# https://www.nature.com/articles/s41467-021-23324-4

rm(list=ls())
sc.sub <- readRDS("Tcells.rds")

library(ProjecTILs)
ref <- load.reference.map() # corresponds to mouse TILs
head(ref[[]])
table(ref$functional.cluster)
table(ref$TILPRED)

library(Seurat)
DefaultAssay(sc.sub) <- "RNA"

# projection:

#projected <- Run.ProjecTILs(sc.sub, ref=ref)
#head(projected[[]])
#plot.projection(ref, projected, linesize = 0.5, pointsize = 0.5)
#plot.statepred.composition(ref, projected, metric = "Percent")
#plot.states.radar(ref, query = projected, min.cells = 30)

# classification:

querydata <- ProjecTILs.classifier(query = sc.sub, ref = ref, split.by = "sample")
head(querydata[[]])
table(querydata$functional.cluster)

p1 <- DimPlot(querydata, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(querydata, group.by = "customclassif", label = F)
p3 <- DimPlot(querydata, group.by = "functional.cluster", label = F)

p1+p2+p3

DimPlot(querydata, group.by = "functional.cluster", label = F, split.by = "sample")


saveRDS(querydata, file = "Tcells.rds")




## cell type proportions ---------------------------------------------------

sc.sub <- readRDS("Tcells.rds")

tab1 <- table(sc.sub$sample, Idents(sc.sub))
prop.table(tab1)
write.csv(tab1, file = "Tcell.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cluster", "Sample","Proportion")
head(tab1b)

ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("")


tab2 <- table(sc.sub$sample, sc.sub$customclassif)
write.csv(tab2, file = "immune.predicted.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("New_clusters", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = New_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Predicted Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2





## markers -----------------------------------------------------------------

sc.sub <- readRDS("Tcells.rds")
DefaultAssay(sc.sub) <- "integrated"

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "markers.Tcell.csv")
write.xlsx(all.markers, file="markers.Tcell.xlsx")








sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample","customclassif","functional.cluster"),
             anno_colors = list(
               alphabet2(20),
               rainbow(14),
               alphabet(14),
               rev(brewer.pal(n = 11, name = "Paired"))),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 4)



# subsetting after re-integration ---------------------------------------------

load(file="cc.int.sct.RData")

DefaultAssay(cc.int.sct) <- "integrated"
head(cc.int.sct[[]])
Idents(cc.int.sct) <- cc.int.sct$sample
levels(Idents(cc.int.sct))

sc.sub <- subset(cc.int.sct, idents = c("DKO_4_Lympho","DKO_8_Lympho","WT_4_Lympho","WT_8_Lympho"))
table(sc.sub$sample)
# DKO_4_Lympho DKO_8_Lympho  WT_4_Lympho  WT_8_Lympho 
# 4465          507         2488         2900 

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 100, dims = 1:30)
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "sample")

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)

Idents(sc.sub)


saveRDS(sc.sub, file = "immune.rds")




# markers -----------------------------------------------------------------

sc.sub <- readRDS("immune.rds")

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "immune.marker.reclustering.csv")
write.xlsx(all.markers, file="immune.marker.reclustering.xlsx")

all.markers <- read.csv(file = "immune.marker.reclustering.csv", row.names = 1)

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample"),
             anno_colors = list(
               alphabet2(24),
               rainbow(4)),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 5)




# automated identification ------------------------------------------------

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])

# ScType
# https://www.nature.com/articles/s41467-022-28803-w
# https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# get cell-type-specific gene sets from our in-built database (DB)
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list = gene_sets_prepare(db_, tissue)
lapply(gs_list, length)
head(gs_list[[1]])
head(gs_list[[2]])


# assign cell types
es.max = sctype_score(scRNAseqData = sc.sub[["integrated"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max[1:10, 1:5]

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(sc.sub@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sc.sub@meta.data[sc.sub@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc.sub@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


# plot

sc.sub@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc.sub@meta.data$customclassif[sc.sub@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sc.sub, reduction = "umap", label = TRUE, label.size = 6, label.box = T, repel = TRUE, group.by = 'customclassif')        
DimPlot(sc.sub, reduction = "umap", label = F, repel = TRUE, group.by = 'customclassif', split.by = "sample")


# bubble plot

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(sc.sub, reduction = "umap", label = TRUE, repel = TRUE, cols = alphabet(23)), gggr, cols = 2)

saveRDS(sc.sub, file = "immune.rds")



# cell type proportions ---------------------------------------------------

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])

tab1 <- table(sc.sub$sample, sc.sub$seurat_clusters)
write.csv(tab1, file = "immune.seurat.recluster.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:24]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cell_type", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Seurat Cluster Proportions") +
  theme_ipsum() +
  xlab("")


tab2 <- table(sc.sub$sample, sc.sub$customclassif)
write.csv(tab2, file = "immune.predicted.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("New_clusters", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = New_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Predicted Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2



# filter out contaminants -------------------------------------------------

rm(list=ls())

sc.sub <- readRDS("immune.rds")
head(sc.sub[[]])
Idents(sc.sub)

DimPlot(sc.sub, label = T, label.size = 6)

FeaturePlot(sc.sub, c("Muc13","Cdh5","Krt14","Krt5","Col4a1","Col6a1","Bmp4"))
VlnPlot(sc.sub, c("Muc13","Cdh5","Krt14","Krt5","Col4a1","Col6a1","Bmp4"))

# cluster 8 : Muc13 expressing cells
# cluster 9 : Stromal cells Bmp4
# cluster 18 : Epi Cdh5
# Cluster 17 : Epi Krt14, Krt5
# Cluster 9, 20, 21, 18 : Stromal Col4a1
# Cluster 9, 20, 21 : Stromal Col6a1

sc.sub <- subset(sc.sub, idents = c('8','9','18','17','20','21'), invert = T)
table(sc.sub$sample)

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 100, dims = 1:30)
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "customclassif", label = T)

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)


saveRDS(sc.sub, file = "immune.2.rds")



## cell type proportions ---------------------------------------------------

sc.sub <- readRDS("immune.2.rds")

tab1 <- table(sc.sub$sample, Idents(sc.sub))
tab1
prop.table(tab1)
write.csv(tab1, file = "immune.ident.cluster.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cluster", "Sample","Proportion")
head(tab1b)

p1 <- ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("")

tab2 <- table(sc.sub$sample, sc.sub$customclassif)
tab2
prop.table(tab2)
write.csv(tab2, file = "immune.ident.cluster.proportions.2.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("Cell_Type", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cell_Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2


## markers -----------------------------------------------------------------

sc.sub <- readRDS("immune.2.rds")
DefaultAssay(sc.sub) <- "integrated"

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 1127 1124 1334 1050 1166 1446 1262 1103 1142 1335 1145 1116 1155 1408 1023  855  943  487  373   
tail(all.markers)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "immune.markers.reclustering.csv")
write.xlsx(all.markers, file="immune.markers.reclustering.xlsx")

all.markers <- read.csv(file = "immune.markers.reclustering.csv", row.names = 1)

sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample","customclassif","functional.cluster"),
             anno_colors = list(
               alphabet2(20),
               rainbow(14),
               alphabet(14),
               rev(brewer.pal(n = 11, name = "Paired"))),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 4)




# T cells -------------------------------------------------

rm(list=ls())

sc.sub <- readRDS("immune.2.rds")
head(sc.sub[[]])
Idents(sc.sub)

Idents(sc.sub) <- sc.sub$customclassif
table(sc.sub$customclassif, sc.sub$sample)

DimPlot(sc.sub, label = T, label.size = 6)

sc.sub <- subset(sc.sub, idents = c('Effector CD4+ T cells'), invert = F)
table(sc.sub$sample)
# DKO_4_Lympho DKO_8_Lympho  WT_4_Lympho  WT_8_Lympho 
# 990           65          264          246 

sc.sub <- SplitObject(sc.sub, split.by = "sample")
sc.sub <- lapply(sc.sub, SCTransform, method = "glmGamPoi")
features  <- SelectIntegrationFeatures(sc.sub, nfeatures = 3000)
sc.sub <- PrepSCTIntegration(sc.sub, anchor.features = features)
sc.sub <- lapply(X = sc.sub, FUN = RunPCA, features = features)
#anchors <- FindIntegrationAnchors(sc.sub, normalization.method = "SCT", anchor.features = features) # this was too slow
anchors <- FindIntegrationAnchors(object.list = sc.sub, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5) # initially used k.anchor = 20, which may be way too high
sc.sub <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 20, dims = 1:30) # reduced k.weight because of low number of cells in one sample
sc.sub <- RunPCA(sc.sub)
sc.sub <- RunUMAP(sc.sub, reduction = "pca", dims = 1:10)
sc.sub <- FindNeighbors(sc.sub, dims = 1:30)
sc.sub <- FindClusters(sc.sub, resolution = 0.5)

p1 <- DimPlot(sc.sub, label = T, label.box = T, label.size = 5)
p2 <- DimPlot(sc.sub, group.by = "customclassif", label = F)

p1+p2

DimPlot(sc.sub, split.by = "sample")#, pt.size = 1)


saveRDS(sc.sub, file = "Tcells.rds")



## ProjecTILs ----

# https://carmonalab.github.io/ProjecTILs_CaseStudies/
# https://github.com/carmonalab/ProjecTILs
# https://www.nature.com/articles/s41467-021-23324-4

rm(list=ls())
sc.sub <- readRDS("Tcells.rds")

library(ProjecTILs)
ref <- load.reference.map() # corresponds to mouse TILs
head(ref[[]])
table(ref$functional.cluster)
table(ref$TILPRED)

library(Seurat)
DefaultAssay(sc.sub) <- "RNA"

# projection:

#projected <- Run.ProjecTILs(sc.sub, ref=ref)
#head(projected[[]])
#plot.projection(ref, projected, linesize = 0.5, pointsize = 0.5)
#plot.statepred.composition(ref, projected, metric = "Percent")
#plot.states.radar(ref, query = projected, min.cells = 30)

# classification:

querydata <- ProjecTILs.classifier(query = sc.sub, ref = ref, split.by = "sample")
head(querydata[[]])
table(querydata$functional.cluster)

p1 <- DimPlot(querydata, label = T, label.box = F, label.size = 5)
p2 <- DimPlot(querydata, group.by = "customclassif", label = F)
p3 <- DimPlot(querydata, group.by = "functional.cluster", label = F)

p1+p2+p3

DimPlot(querydata, group.by = "functional.cluster", label = F, split.by = "sample")


saveRDS(querydata, file = "Tcells.rds")




## cell type proportions ---------------------------------------------------

sc.sub <- readRDS("Tcells.rds")

tab1 <- table(sc.sub$sample, Idents(sc.sub))
prop.table(tab1)
write.csv(tab1, file = "Tcell.proportions.csv")

colors <- as.vector(glasbey(n=24))
cluster.colors = colors[1:20]
condition.colors = c("brown","orange","yellow","pink")

tab1b <- as.data.frame(t(tab1))
colnames(tab1b) <- c("Cluster", "Sample","Proportion")
head(tab1b)

ggplot(tab1b, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Cluster Proportions") +
  theme_ipsum() +
  xlab("")


tab2 <- table(sc.sub$sample, sc.sub$customclassif)
write.csv(tab2, file = "immune.predicted.celltype.proportions.csv")

tab2b <- as.data.frame(t(tab2))
colnames(tab2b) <- c("New_clusters", "Sample","Proportion")
head(tab2b)

p2 <- ggplot(tab2b, aes(x = Sample, y = Proportion, fill = New_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster.colors) +
  ggtitle("Predicted Cell Type Proportions") +
  theme_ipsum() +
  xlab("")

p1+p2





## markers -----------------------------------------------------------------

sc.sub <- readRDS("Tcells.rds")
DefaultAssay(sc.sub) <- "integrated"

all.markers <- FindAllMarkers(sc.sub, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
table(all.markers$cluster)

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(all.markers, file = "markers.Tcell.csv")
write.xlsx(all.markers, file="markers.Tcell.xlsx")








sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 5, wt = avg_log2FC)

plot_heatmap(dataset = sc.sub,  
             markers = unique(sel.markers$gene),
             anno_var = c("seurat_clusters","sample","customclassif","functional.cluster"),
             anno_colors = list(
               alphabet2(20),
               rainbow(14),
               alphabet(14),
               rev(brewer.pal(n = 11, name = "Paired"))),
             row_font_size = 10
)


sel.markers <- group_by(all.markers, cluster) %>% top_n(n = 1, wt = avg_log2FC)
sel.markers <- unique(sel.markers$gene)

FeaturePlot(sc.sub, features = sel.markers, ncol = 4)




# end ---------------------------------------------------------------------







