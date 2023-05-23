
# scRNAseq ----

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell/results/scRNAseq/cellchat/")

## role of CD137 ----

check <- read.csv("all.celltypes/C1D1_interactions.csv", row.names = 1)
check <- read.csv("all.celltypes/C3D1_interactions.csv", row.names = 1)
head(check)

check[check$ligand=="TNFSF9", ]
check[grep("CD137", check$pathway_name), ]


## nature of increased MHC1 interactions ----

check <- read.csv("Tcells/C1D1_interactions.csv", row.names = 1)
check <- read.csv("Tcells/C3D1_interactions.csv", row.names = 1)
head(check)

check[check$ligand=="TNFSF9", ]
mhc1 <- check[grep("MHC-I", check$pathway_name), ]
dim(mhc1)
head(mhc1)

# 44 MHC1 interactions for C1D1 and 133 for C3D1


## Monocyte:DC switch ----

check1 <- read.csv("all.celltypes/C1D1_interactions.csv", row.names = 1)
check1 <- read.csv("all.celltypes/C3D1_interactions.csv", row.names = 1)
head(check1)

table(check1$source)
table(check1$target)

# 190/295 monocyte interactions before therapy (vs.0 after)
# 114/107 DC interactions after therapy (vs.0 before)

mono <- check1[check1$source=="Monocytes", ]
mono <- check1[check1$target=="Monocytes", ]
head(mono)
table(mono$ligand)
table(mono$receptor)

dc <- check1[check1$source=="Dendritic cells", ]
dc <- check1[check1$target=="Dendritic cells", ]
head(dc)
table(dc$ligand)
table(dc$receptor)


## role of endothelial cells ----

check1 <- read.csv("all.celltypes/C1D1_interactions.csv", row.names = 1)
check2 <- read.csv("all.celltypes/C3D1_interactions.csv", row.names = 1)
head(check1)

table(check1$source) # 481 endothelial
table(check1$target) # 337
table(check2$source) # 603 endothelial
table(check2$target) # 456

endo.source1 <- check1[check1$source=="Endothelial cells", ]
endo.target1 <- check1[check1$target=="Endothelial cells", ]
endo.source2 <- check2[check2$source=="Endothelial cells", ]
endo.target2 <- check2[check2$target=="Endothelial cells", ]

setdiff(endo.source1$interaction_name, endo.source2$interaction_name)
# "PDGFD_PDGFRB" "APLN_APLNR"   "CXCL2_ACKR1"  "CXCL2_CXCR2"  "SELE_CEACAM1" "EFNA1_EPHA2"
setdiff(endo.source2$interaction_name, endo.source1$interaction_name)
#"PDGFB_PDGFRA"       "VEGFC_VEGFR3"       "VEGFC_VEGFR2"       "PGF_VEGFR1"        
#[5] "VEGFC_VEGFR2R3"     "CCL14_CCR1"         "CCL14_ACKR1"        "CXCL8_CXCR1"       
#[9] "MDK_SDC1"           "MDK_SDC2"           "MDK_SDC4"           "MDK_ITGA4_ITGB1"   
#[13] "MDK_ITGA6_ITGB1"    "MDK_LRP1"           "MDK_NCL"            "POSTN_ITGAV_ITGB5" 
#[17] "EDN1_EDNRA"         "EDN1_EDNRB"         "FN1_ITGA3_ITGB1"    "FN1_ITGA8_ITGB1"   
#[21] "LAMA4_ITGA3_ITGB1"  "LAMA5_ITGA3_ITGB1"  "LAMB1_ITGA3_ITGB1"  "LAMB2_ITGA3_ITGB1" 
#[25] "LAMC1_ITGA3_ITGB1"  "THBS1_ITGA3_ITGB1"  "COL4A1_ITGA3_ITGB1" "COL4A2_ITGA3_ITGB1"
#[29] "COL6A1_ITGA3_ITGB1" "COL6A2_ITGA3_ITGB1" "TNXB_ITGA8_ITGB1"   "TNXB_ITGA9_ITGB1"  
#[33] "COL4A1_ITGA9_ITGB1" "COL4A2_ITGA9_ITGB1" "COL6A1_ITGA9_ITGB1" "COL6A2_ITGA9_ITGB1"
#[37] "LAMA4_ITGA9_ITGB1"  "LAMA5_ITGA9_ITGB1"  "LAMB1_ITGA9_ITGB1"  "LAMB2_ITGA9_ITGB1" 
#[41] "LAMC1_ITGA9_ITGB1"  "TNXB_SDC1"          "TNXB_SDC4"          "AGRN_DAG1"         
#[45] "EFNA1_EPHA3"        "EFNA1_EPHA4"        "EFNB1_EPHA4"        "EFNB1_EPHB4"       
#[49] "EFNB2_EPHA4"        "JAM2_ITGA3_ITGB1"   "F11R_F11R"          "F11R_JAM2"         
#[53] "F11R_JAM3"          "HLA-A_CD8B"         "HLA-B_CD8B"         "HLA-C_CD8B"        
#[57] "HLA-E_CD8B"         "HLA-F_CD8B"         "HLA-E_KLRC2"        "HLA-F_LILRB1"      
#[61] "HLA-E_CD94:NKG2C"   "HLA-DPB1_CD4"       "HLA-DMA_CD4"        "HLA-DQB1_CD4"      
#[65] "PVR_TIGIT"          "SEMA4C_PLXNB2"      "SEMA6A_PLXNA2"      "ITGA9_ITGB1_VCAM1" 


pathways1 <- unique(c(endo.source1$pathway_name, endo.target1$pathway_name)) # 60
pathways2 <- unique(c(endo.source2$pathway_name, endo.target2$pathway_name)) # 66

setdiff(pathways1, pathways2)
# "APELIN" "CEACAM" "IGF"    "SPP1" 
setdiff(pathways2, pathways1)
# "PERIOSTIN" "EDN"       "TENASCIN"  "AGRN"      "PVR"       "ncWNT"     "OSM"      
#[8] "LIGHT"     "CD96"      "SEMA7"

endo.source2[endo.source2$pathway_name=="PERIOSTIN", ]
endo.target2[endo.target2$pathway_name=="PERIOSTIN", ]

## endothelial cells post-therapy target tumor EMT high cells through POSTN:ITGAV_ITGB5 interaction

endo.source2[endo.source2$pathway_name=="CD96", ]
endo.target2[endo.target2$pathway_name=="CD96", ]

degs <- read.csv("/home/hh/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell/results/scRNAseq/seurat/all.cells/DEGs.Endothelial cells.csv", row.names = 1)
head(degs)
rownames(degs)
table(degs$avg_log2FC>0)


## interactions with tumor cells ----

check1 <- read.csv("all.celltypes/C1D1_interactions.csv", row.names = 1)
check2 <- read.csv("all.celltypes/C3D1_interactions.csv", row.names = 1)
head(check1)

table(check1$target)
table(check2$target)

check1 <- check1[check1$target=="Tumor_EMT.high", ]
check2 <- check2[check2$target=="Tumor_EMT.high", ]

setdiff(check1$pathway_name, check2$pathway_name)
# "IGF"    "IL6"    "SPP1"   "CEACAM" "EPHA"
setdiff(check2$pathway_name, check1$pathway_name)
# "ncWNT"     "LIFR"      "OSM"       "LIGHT"     "PERIOSTIN" "TENASCIN"  "AGRN" 

unique.paths <- setdiff(check2$pathway_name, check1$pathway_name)
sel.check2 <- check2[check2$pathway_name%in%unique.paths, ]


## CAFs ----

degs <- read.csv("/home/hh/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell/results/scRNAseq/seurat/all.cells/DEGs.CAFs.csv", row.names = 1)
head(degs)
rownames(degs)
table(degs$avg_log2FC>0)


# spatial ----

setwd("~/Dropbox/BioInfo/Colabs/Mehlen/NP137_single.cell/results/spatial/cellchat/")

## endothelial ----

# as target patient 1

check1 <- read.csv("interactions_P034_Pre.csv", row.names = 1)
check2 <- read.csv("interactions_P034_Post.csv", row.names = 1)
head(check1)

table(check1$target) # 46 endothelial
table(check2$target) # 683

check1 <- check1[check1$target=="Endothelial", ]
check2 <- check2[check2$target=="Endothelial", ]

setdiff(check1$pathway_name, check2$pathway_name)
# "OSM"   "TWEAK" "MK"    "PTN"   "NT"    "APP"   "OCLN" 
setdiff(check2$pathway_name, check1$pathway_name)
#[1] "BMP"        "ACTIVIN"    "WNT"        "ncWNT"      "EGF"        "NRG"        "FGF"        "APELIN"     "HH"         "CCL"        "CXCL"       "MIF"        "IL1"        "CSF"        "IL16"       "TNF"        "VEGI"      
#[18] "VISFATIN"   "ANGPTL"     "PERIOSTIN"  "BRADYKININ" "COMPLEMENT" "PARs"       "HGF"        "SEMA3"      "CALCR"      "GAS"        "PROS"       "CHEMERIN"   "RELN"       "VTN"        "HSPG"       "CD45"       "CD46"      
#[35] "CDH"        "CDH5"       "EPHA"       "EPHB"       "ESAM"       "JAM"        "MHC-I"      "MHC-II"     "SELE"       "SELPLG"     "SEMA4"      "SEMA5" 

unique.paths <- setdiff(check2$pathway_name, check1$pathway_name)
sel.check2 <- check2[check2$pathway_name%in%unique.paths, ]
table(sel.check2$source)
table(sel.check2$target)


# as source patient 1

check1 <- read.csv("interactions_P034_Pre.csv", row.names = 1)
check2 <- read.csv("interactions_P034_Post.csv", row.names = 1)

table(check1$source) # 46
table(check2$source) # 706

check1 <- check1[check1$source=="Endothelial", ]
check2 <- check2[check2$source=="Endothelial", ]

setdiff(check1$pathway_name, check2$pathway_name)
#  "OSM"   "TWEAK" "PTN"   "NT"    "APP"   "OCLN"  "SEMA6" 
setdiff(check2$pathway_name, check1$pathway_name)
#[1] "BMP"        "ACTIVIN"    "WNT"        "ncWNT"      "FGF"        "APELIN"     "CCL"        "CXCL"       "MIF"        "CSF"        "TNF"        "VEGI"       "VISFATIN"   "ANGPTL"     "PERIOSTIN"  "BRADYKININ" "COMPLEMENT"
#[18] "PARs"       "HGF"        "SEMA3"      "GAS"        "PROS"       "CHEMERIN"   "VTN"        "HSPG"       "CD45"       "CD46"       "CDH"        "CDH5"       "CNTN"       "EPHA"       "EPHB"       "ESAM"       "JAM"       
#[35] "MHC-I"      "MHC-II"     "SELE"       "SELPLG"     "SEMA5"      "SN"         "VCAM" 

unique.paths <- setdiff(check2$pathway_name, check1$pathway_name)
sel.check2 <- check2[check2$pathway_name%in%unique.paths, ]
table(sel.check2$source)
table(sel.check2$target)



# as target patient 2

check1 <- read.csv("interactions_P039_Pre.csv", row.names = 1)
check2 <- read.csv("interactions_P039_Post.csv", row.names = 1)
head(check1)

table(check1$target) # 524
table(check2$target) # 1642
table(check1$source) # 454
table(check2$source) # 1429

check1 <- check1[check1$target=="Endothelial", ]
check2 <- check2[check2$target=="Endothelial", ]

setdiff(check1$pathway_name, check2$pathway_name)
# "IL6"       "TNF"       "TWEAK"     "RELN"      "CD34"      "CD46"      "CD6"       "DESMOSOME" "NECTIN"    "OCLN" 
setdiff(check2$pathway_name, check1$pathway_name)
#[1] "NRG"        "IGF"        "APELIN"     "IL2"        "VEGI"       "EDA"        "COMPLEMENT" "TAC"        "APJ"        "PROS"       "PTH"        "NPNT"       "ADGRE5"     "APP"        "CD86"       "CD96"       "CDH"       
#[18] "CDH1"       "ICAM"       "JAM"        "NCAM"       "NEGR"       "THY1" 

unique.paths <- setdiff(check2$pathway_name, check1$pathway_name)
sel.check2 <- check2[check2$pathway_name%in%unique.paths, ]
table(sel.check2$source)
table(sel.check2$target)

sel.check2[sel.check2$source == "Tumor cells", ]


# as source patient 2

check1 <- read.csv("interactions_P039_Pre.csv", row.names = 1)
check2 <- read.csv("interactions_P039_Post.csv", row.names = 1)

check1 <- check1[check1$source=="Endothelial", ]
check2 <- check2[check2$source=="Endothelial", ]

setdiff(check1$pathway_name, check2$pathway_name)
#   "EGF"       "IL6"       "TWEAK"     "PARs"      "AGRN"      "HSPG"      "ALCAM"     "CD34"      "DESMOSOME" "NECTIN"    "OCLN"      "SELL" 
setdiff(check2$pathway_name, check1$pathway_name)
#"GDF"    "APELIN" "IL2"    "EDA"    "PTN"    "NMU"    "SEMA3"  "CALCR"  "APJ"    "PSAP"   "ADGRE5" "APP"    "CDH"    "CDH1"   "CNTN"   "CSPG4"  "JAM"    "L1CAM"  "NCAM"   "NEGR"   "SELPLG" "SEMA4"  "SEMA6"  "VISTA" 

unique.paths <- setdiff(check2$pathway_name, check1$pathway_name)
sel.check2 <- check2[check2$pathway_name%in%unique.paths, ]
table(sel.check2$source)
table(sel.check2$target)

sel.check2[sel.check2$target == "Tumor cells", ]




# end ----
sessionInfo()
