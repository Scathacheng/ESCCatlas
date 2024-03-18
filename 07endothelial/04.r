##R 
##single cell data of ESCC
##step4: cell-cell communication analysis
##Jiacheng Dai #daicy0424@gmail.com
##2022/08/18

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(CellChat) #v1.0.0
library(ggalluvial)
library(ggplot2)


#load data
setwd("01ECatlas/04CCI")
data = readRDS("../02EC_cluster/02ECdataset.sct.rds")
meta = read.csv("../02EC_cluster/02ECdataset_Meta.csv",header=T,row.names=1)
meta = meta[meta$Doublet=="Singlet",]

meta.tumor = meta[meta$type=="Tumor",]
meta.adj = meta[meta$type=="Peri",]
meta.dist = meta[meta$type=="Normal",]
data.matrix.tumor = data@assays$SCT@data[,rownames(meta.tumor)]
data.matrix.adj = data@assays$SCT@data[,rownames(meta.adj)]
data.matrix.dist = data@assays$SCT@data[,rownames(meta.dist)]

###Cellchat
data.input = data.matrix.tumor
meta = meta.tumor
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

#add cell information
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
#set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)  # is a plot
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use Secreted Signaling for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
#preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
### inference of cell-cell communication network
#compute the communication probability and infer cellular commincation network
cellchat <- computeCommunProb(cellchat,raw.use=TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


#选择cell-cell communication 
df.net <- subsetCommunication(cellchat, sources.use = c(2,3,4,5,6,7,8,9,10), targets.use = c(1,11,12,13,14))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(2,3,4,5,6,7,8,9,10), targets.use = c(1,11,12,13,14), remove.isolate=FALSE)
### visualization of cell-cell communication network
pdf("bubble_dist.pdf",h=18,w=20)
netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6,7,8,9,10), targets.use = c(1,11,12,13,14), remove.isolate = FALSE)
dev.off()
pdf("chord_dist.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(2,3,4,5,6,7,8,9,10), targets.use = c(1,11,12,13,14), lab.cex = 0.5,legend.pos.y = 30)
dev.off()


#selected clusters
data = readRDS("../02EC_cluster/02ECdataset.sct.rds")
meta = read.csv("../02EC_cluster/02ECdataset_Meta.csv",header=T,row.names=1)
meta = meta[meta$Doublet=="Singlet",]

list_cluster1 = c(15,17,18,20,23,0,4,14,8,9,3)
list_cluster2 = c("KRT+ unknown-C15", "fibroblast-like-C17", "immune-like-C18", "SMC-like-C20", "Macro-like-C23","HLA+ vein-C0", "CCL2+ vein-C4", "HSP+ vein-C14", "FABP4+ profEC-C8", "SOCS3+ profEC-C9","COL4A1+ TEC-C3")

meta1 = meta[meta$seurat_clusters_EC == list_cluster1[1],]
meta1$labels = list_cluster2[1]
meta2 = meta[meta$seurat_clusters_EC == list_cluster1[2],]
meta2$labels = list_cluster2[2]
meta3 = meta[meta$seurat_clusters_EC == list_cluster1[3],]
meta3$labels = list_cluster2[3]
meta4 = meta[meta$seurat_clusters_EC == list_cluster1[4],]
meta4$labels = list_cluster2[4]
meta5 = meta[meta$seurat_clusters_EC == list_cluster1[5],]
meta5$labels = list_cluster2[5]
meta6 = meta[meta$seurat_clusters_EC == list_cluster1[6],]
meta6$labels = list_cluster2[6]
meta7 = meta[meta$seurat_clusters_EC == list_cluster1[7],]
meta7$labels = list_cluster2[7]
meta8 = meta[meta$seurat_clusters_EC == list_cluster1[8],]
meta8$labels = list_cluster2[8]
meta9 = meta[meta$seurat_clusters_EC == list_cluster1[9],]
meta9$labels = list_cluster2[9]
meta10 = meta[meta$seurat_clusters_EC == list_cluster1[10],]
meta10$labels = list_cluster2[10]
meta11 = meta[meta$seurat_clusters_EC == list_cluster1[11],]
meta11$labels = list_cluster2[11]

meta_new = rbind(meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10, meta11)
meta.tumor<-meta_new[which(meta_new$type=="Tumor"),]
meta.peri<-meta_new[which(meta_new$type=="Peri"),]
meta.normal<-meta_new[which(meta_new$type=="Normal"),]
data.matrix<-data@assays$SCT@data
data.matrix.tumor<-data.matrix[,rownames(meta.tumor)]
data.matrix.peri<-data.matrix[,rownames(meta.peri)]
data.matrix.normal<-data.matrix[,rownames(meta.normal)]

data.input = data.matrix.normal
meta = meta.normal
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

#选择cell-cell communication 
df.net <- subsetCommunication(cellchat, sources.use = c(4,7,8,9,10), targets.use = c(1,2,3,5,6,11))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(4,7,8,9,10), targets.use = c(1,2,3,5,6,11), remove.isolate=FALSE)
### visualization of cell-cell communication network
pdf("bubble_dist_select.pdf",h=8,w=12)
netVisual_bubble(cellchat, sources.use = c(4,7,8,9,10), targets.use = c(1,2,3,5,6,11), remove.isolate = FALSE)
dev.off()
pdf("chord_dist_select.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(4,7,8,9,10), targets.use = c(1,2,3,5,6,11), lab.cex = 0.5,legend.pos.y = 30)
dev.off()


# 看其他T细胞等细胞类型的 03test
options(stringsAsFactors=TRUE)
meta = read.csv("01ESCCdataset_Meta.csv",header=T,row.names=1)
data = readRDS("01ESCCdataset.sct.rds")
clusters.T = meta[meta$seurat_clusters %in% c(0,2,21,23,24,26,28),]
clusters.T$celltype = "T cell"
clusters.M = meta[meta$seurat_clusters %in% c(13,30),]
clusters.M$celltype = "Mast cell"
clusters.B = meta[meta$seurat_clusters %in% c(1,32,33,39,40,41,42,44,45,46,49),]
clusters.B$celltype = "B cell"
clusters.epi = meta[meta$seurat_clusters %in% c(8,14,18,19,25,35,48),]
clusters.epi$celltype = "Epithelial cell"
clusters.myeloid = meta[meta$seurat_clusters %in% c(12,16,17,34,36,47),]
clusters.myeloid$celltype = "Myeloid cell"
clusters.E = meta[meta$seurat_clusters %in% c(3,9,50,51),]
clusters.E$celltype = "Endothelial cell"
clusters.L = meta[meta$seurat_clusters %in% c(31),]
clusters.L$celltype = "Lymphatic Endothelial cell"
clusters.F = meta[meta$seurat_clusters %in% c(4,5,6,10,11,29,38),]
clusters.F$celltype = "Fibroblast"
clusters.smc = meta[meta$seurat_clusters %in% c(7,15,20,22),]
clusters.smc$celltype = "Smooth Muscle Cell"
clusters.o = meta[meta$seurat_clusters %in% c(37,43),]
clusters.o$celltype = "other"
clusters = rbind(clusters.T, clusters.M, clusters.B, clusters.epi, clusters.myeloid, clusters.E, clusters.L, clusters.F, clusters.smc, clusters.o)

meta.tumor<-clusters[which(clusters$type=="Tumor"),]
meta.peri<-clusters[which(clusters$type=="Peri"),]
meta.normal<-clusters[which(clusters$type=="Normal"),]
data.matrix<-data@assays$SCT@data
data.matrix.tumor<-data.matrix[,rownames(meta.tumor)]
data.matrix.peri<-data.matrix[,rownames(meta.peri)]
data.matrix.normal<-data.matrix[,rownames(meta.normal)]

data.input = data.matrix.normal
meta = meta.normal
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

#选择cell-cell communication 
df.net <- subsetCommunication(cellchat, sources.use = c(1,3,4,5,6,7,8,9,10), targets.use = c(2))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(1,3,4,5,6,7,8,9,10), targets.use = c(2), remove.isolate=FALSE)
### visualization of cell-cell communication network
pdf("03test/bubble_dist.pdf",h=8,w=12)
netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8,9,10), targets.use = c(2), remove.isolate = FALSE)
dev.off()
pdf("03test/chord_dist.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(1,3,4,5,6,7,8,9,10), targets.use = c(2), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

# 04test
meta1 = read.csv("../01total_cluster/01ESCCdataset_Meta.csv",header=T,row.names=1)
clusters.T = meta1[meta1$seurat_clusters %in% c(0,2,21,23,24,26,28),]
clusters.T$celltype = "T cell"
clusters.F = meta1[meta1$seurat_clusters %in% c(4,5,6,10,11,29,38),]
clusters.F$celltype = "Fibroblast"

meta2 = read.csv("../02EC_cluster/02ECdataset_Meta.csv",header=T,row.names=1)
meta = meta2[meta2$Doublet=="Singlet",]
list_cluster1 = c(0,2,4,12,14,16,8,9,3,19)
list_cluster2 = c("HLA+ vein-C0","RGCC+ EC-C2","CCL2+ vein-C4","ISG15+ vein-C12","HSP+ vein-C14","ACTG1+ vein-C16","FABP4+ EC-C8", "SOCS3+ EC-C9","COL4A1+ TEC-C3","APLNR+ TEC-C19")
meta.1 = meta[meta$seurat_clusters_EC == list_cluster1[1],]
meta.1$celltype = list_cluster2[1]
meta.2 = meta[meta$seurat_clusters_EC == list_cluster1[2],]
meta.2$celltype = list_cluster2[2]
meta.3 = meta[meta$seurat_clusters_EC == list_cluster1[3],]
meta.3$celltype = list_cluster2[3]
meta.4 = meta[meta$seurat_clusters_EC == list_cluster1[4],]
meta.4$celltype = list_cluster2[4]
meta.5 = meta[meta$seurat_clusters_EC == list_cluster1[5],]
meta.5$celltype = list_cluster2[5]
meta.6 = meta[meta$seurat_clusters_EC == list_cluster1[6],]
meta.6$celltype = list_cluster2[6]
meta.7 = meta[meta$seurat_clusters_EC == list_cluster1[7],]
meta.7$celltype = list_cluster2[7]
meta.8 = meta[meta$seurat_clusters_EC == list_cluster1[8],]
meta.8$celltype = list_cluster2[8]
meta.9 = meta[meta$seurat_clusters_EC == list_cluster1[9],]
meta.9$celltype = list_cluster2[9]
meta.10 = meta[meta$seurat_clusters_EC == list_cluster1[10],]
meta.10$celltype = list_cluster2[10]

meta_new = rbind(meta.1, meta.2, meta.3, meta.4, meta.5, meta.6, meta.7, meta.8, meta.9, meta.10)
meta_new = meta_new[,-12:-13]
meta_total = rbind(clusters.T, clusters.F, meta_new)

meta.tumor<-meta_total[which(meta_total$type=="Tumor"),]
meta.peri<-meta_total[which(meta_total$type=="Peri"),]
meta.normal<-meta_total[which(meta_total$type=="Normal"),]
data = readRDS("../01total_cluster/01ESCCdataset.sct.rds")
data.matrix<-data@assays$SCT@data
data.matrix.tumor<-data.matrix[,rownames(meta.tumor)]
data.matrix.peri<-data.matrix[,rownames(meta.peri)]
data.matrix.normal<-data.matrix[,rownames(meta.normal)]

data.input = data.matrix.normal
meta = meta.normal
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

#选择cell-cell communication 
df.net <- subsetCommunication(cellchat, sources.use = c(6,12), targets.use = c(1,2,3,4,5,7,8,9,10,11))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(6,12), targets.use = c(1,2,3,4,5,7,8,9,10,11), remove.isolate=FALSE)
### visualization of cell-cell communication network
pdf("04test/bubble_tumor.pdf",h=8,w=12)
netVisual_bubble(cellchat, sources.use = c(6,12), targets.use = c(1,2,3,4,5,7,8,9,10,11), remove.isolate = FALSE)
dev.off()
pdf("04test/chord_tumor.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(6,12), targets.use = c(1,2,3,4,5,7,8,9,10,11), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

#05 test
#! /usr/bin/env Rscript

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(CellChat) #v1.0.0
library(ggalluvial)
library(ggplot2)

meta1 = read.csv("../01total_cluster/01ESCCdataset_Meta.csv",header=T,row.names=1)
clusters.T = meta1[meta1$seurat_clusters %in% c(0,2,21,23,24,26,28),]
clusters.T$celltype = "T cell"
clusters.F = meta1[meta1$seurat_clusters %in% c(4,5,6,10,11,29,38),]
clusters.F$celltype = "Fibroblast"
clusters.epi = meta1[meta1$seurat_clusters %in% c(8,14,18,19,25,35,48),]
clusters.epi$celltype = "Epithelial cell"
clusters.myeloid = meta1[meta1$seurat_clusters %in% c(12,16,17,34,36,47),]
clusters.myeloid$celltype = "Myeloid cell"
clusters.B = meta1[meta1$seurat_clusters %in% c(1,32,33,39,40,41,42,44,45,46,49),]
clusters.B$celltype = "B cell"

meta2 = read.csv("../02EC_cluster/02ECdataset_Meta.csv",header=T,row.names=1)
meta = meta2[meta2$Doublet=="Singlet",]
list_cluster1 = c(0,2,4,10,12,14,16,22,8,9,3,19,24,5,6,13,1)
list_cluster2 = c("HLA+ vein","RGCC+ EC","CCL2+ vein","CPE+ vein","ISG15+ vein","HSP+ EC","ACTG1+ vein","SELE+ vein","FABP4+ EC", "SOCS3+ EC","COL4A1+ TEC","APLNR+ TEC","TOP2A+ TEC","JAG1+ artery","CAV1+ cap","MGP+ artery","Qsc EC")
meta.1 = meta[meta$seurat_clusters_EC == list_cluster1[1],]
meta.1$celltype = list_cluster2[1]
meta.2 = meta[meta$seurat_clusters_EC == list_cluster1[2],]
meta.2$celltype = list_cluster2[2]
meta.3 = meta[meta$seurat_clusters_EC == list_cluster1[3],]
meta.3$celltype = list_cluster2[3]
meta.4 = meta[meta$seurat_clusters_EC == list_cluster1[4],]
meta.4$celltype = list_cluster2[4]
meta.5 = meta[meta$seurat_clusters_EC == list_cluster1[5],]
meta.5$celltype = list_cluster2[5]
meta.6 = meta[meta$seurat_clusters_EC == list_cluster1[6],]
meta.6$celltype = list_cluster2[6]
meta.7 = meta[meta$seurat_clusters_EC == list_cluster1[7],]
meta.7$celltype = list_cluster2[7]
meta.8 = meta[meta$seurat_clusters_EC == list_cluster1[8],]
meta.8$celltype = list_cluster2[8]
meta.9 = meta[meta$seurat_clusters_EC == list_cluster1[9],]
meta.9$celltype = list_cluster2[9]
meta.10 = meta[meta$seurat_clusters_EC == list_cluster1[10],]
meta.10$celltype = list_cluster2[10]
meta.11 = meta[meta$seurat_clusters_EC == list_cluster1[11],]
meta.11$celltype = list_cluster2[11]
meta.12 = meta[meta$seurat_clusters_EC == list_cluster1[12],]
meta.12$celltype = list_cluster2[12]
meta.13 = meta[meta$seurat_clusters_EC == list_cluster1[13],]
meta.13$celltype = list_cluster2[13]
meta.14 = meta[meta$seurat_clusters_EC == list_cluster1[14],]
meta.14$celltype = list_cluster2[14]
meta.15 = meta[meta$seurat_clusters_EC == list_cluster1[15],]
meta.15$celltype = list_cluster2[15]
meta.16 = meta[meta$seurat_clusters_EC == list_cluster1[16],]
meta.16$celltype = list_cluster2[16]
meta.17 = meta[meta$seurat_clusters_EC == list_cluster1[17],]
meta.17$celltype = list_cluster2[17]

meta_new = rbind(meta.1, meta.2, meta.3, meta.4, meta.5, meta.6, meta.7, meta.8, meta.9, meta.10, meta.11, meta.12, meta.13, meta.14, meta.15, meta.16, meta.17)
meta_new = meta_new[,-12:-13]
meta_total = rbind(clusters.T, clusters.F, clusters.epi, clusters.B, clusters.myeloid, meta_new)

meta.tumor<-meta_total[which(meta_total$type=="Tumor"),]
meta.peri<-meta_total[which(meta_total$type=="Peri"),]
meta.normal<-meta_total[which(meta_total$type=="Normal"),]
data = readRDS("../01total_cluster/01ESCCdataset.sct.rds")
data.matrix<-data@assays$SCT@data
data.matrix.tumor<-data.matrix[,rownames(meta.tumor)]
data.matrix.peri<-data.matrix[,rownames(meta.peri)]
data.matrix.normal<-data.matrix[,rownames(meta.normal)]

data.input = data.matrix.normal
meta = meta.normal
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

#add cell information
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
#set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
#showDatabaseCategory(CellChatDB)  # is a plot
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use Secreted Signaling for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
#preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
### inference of cell-cell communication network
#compute the communication probability and infer cellular commincation network
cellchat <- computeCommunProb(cellchat,raw.use=TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#选择cell-cell communication 
df.net <- subsetCommunication(cellchat, sources.use = c(3,8,10,16,21), targets.use = c(1,2,4,5,6,7,9,11,12,13,14,15,17,18,19,20,22))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(3,8,10,16,21), targets.use = c(1,2,4,5,6,7,9,11,12,13,14,15,17,18,19,20,22), remove.isolate=FALSE)
### visualization of cell-cell communication network
pdf("07test/bubble_dist.pdf",h=8,w=20)
netVisual_bubble(cellchat, sources.use = c(3,8,10,16,21), targets.use = c(1,2,4,5,6,7,9,11,12,13,14,15,17,18,19,20,22), remove.isolate = FALSE)
dev.off()
pdf("07test/chord_dist.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(3,8,10,16,21), targets.use = c(1,2,4,5,6,7,9,11,12,13,14,15,17,18,19,20,22), lab.cex = 0.5,legend.pos.y = 30)
dev.off()
saveRDS(cellchat, "07test/cellchat_dist.rds")

cellchat = netAnalysis_computeCentrality(cellchat, slot.name="netP") 

pathways.show=c("VEGF")
pdf("07test/netvisual_agg_VEGF_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_VEGF_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_adj_VEGF.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("PDGF")
pdf("07test/netvisual_agg_PDGF_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_PDGF_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_PDGF_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("CCL")
pdf("07test/netvisual_agg_CCL_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_CCL_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_CCL_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("FGF")
pdf("07test/netvisual_agg_FGF_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_FGF_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_FGF_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("TGFb")
pdf("07test/netvisual_agg_TGFb_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_TGFb_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_TGFb_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("BMP")
pdf("07test/netvisual_agg_BMP_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_BMP_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_BMP_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("CXCL")
pdf("07test/netvisual_agg_CXCL_adj.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("07test/netvisual_hmap_CXCL_adj.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("07test/sgnRole_CXCL_adj.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()





