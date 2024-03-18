#R
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(CellChat) #v1.0.0
library(ggalluvial)
library(ggplot2)

#load EC data
data1 = readRDS("../02EC_cluster/new0423/02ECdataset.16580.rds")
#load fibroblast data
data2 = readRDS("../08fibroblast/08FB.sct.rds")
#load SMC data
data3 = readRDS("../08smc/08SMC.sct.rds") 

genelist = Reduce(intersect, list(rownames(data1),rownames(data2),rownames(data3)))
#17980

#1 TEC interaction
meta1 = data1@meta.data[data1@meta.data$celltype %in% c("Artery_UNC5B+","tip cell_COL4A1+","stalk cell_RGCC+","Vein_SELE+"),c("celltype","seurat_clusters")]
meta2 = data2@meta.data[data2@meta.data$celltype %in% c("CAF_1","CAF_2","CAF_3"),c("celltype","seurat_clusters")]
meta3 = data3@meta.data[data3@meta.data$celltype %in% c("T-Pericyte","T-SMC"),c("celltype","seurat_clusters")]

#2 immune interaction
meta1 = data1@meta.data[data1@meta.data$celltype %in% c("Vein_CPE+","Vein_HLA+","Vein_IL6+","Vein_SELE+","Vein_ISG15+","PCV_ACTG1+","EC_HSP+"),c("celltype","seurat_clusters","tissue")]
meta2 = data2@meta.data[data2@meta.data$celltype %in% c("FB_1","FB_2","FB_3","FB_4","FB_5","FB_6","FB_7"),c("celltype","seurat_clusters","tissue")]
meta3 = data3@meta.data[data3@meta.data$celltype %in% c("SMC_1","SMC_2","SMC_3","Pericyte_1","Pericyte_2","Pericyte_3","Pericyte_4","Pericyte_5","Pericyte_6"),c("celltype","seurat_clusters","tissue")]

matrix1 = data1@assays$SCT@data[genelist,rownames(meta1)]
matrix2 = data2@assays$SCT@data[genelist,rownames(meta2)]
matrix3 = data3@assays$SCT@data[genelist,rownames(meta3)]

###Cellchat
data.input = cbind(matrix1, matrix2, matrix3)
meta = rbind(meta1, meta2, meta3)

Distant = meta[meta$tissue=="Distant",]
data.Distant = data.input[,rownames(Distant)]
cellchat <- createCellChat(object = data.Distant, meta = Distant, group.by = "celltype")

#add cell information
cellchat <- addMeta(cellchat, meta = Distant)
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
df.net <- subsetCommunication(cellchat, sources.use = c(2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18), targets.use = c(1,9,19,20,21,22,23))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18), targets.use = c(1,9,19,20,21,22,23), remove.isolate=FALSE)

### visualization of cell-cell communication network
pdf("bubble_dist.pdf",h=12,w=15)
netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18), targets.use = c(1,9,19,20,21,22,23), remove.isolate = FALSE)
dev.off()
pdf("chord_dist.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18), targets.use = c(1,9,19,20,21,22,23), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

#highlight VEGF pathway
cellchat = netAnalysis_computeCentrality(cellchat, slot.name="netP") 

pathways.show=c("VEGF")
pdf("netvisual_agg_VEGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("netvisual_hmap_VEGF.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("sgnRole_VEGF.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("PDGF")
pdf("netvisual_agg_PDGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("netvisual_hmap_PDGF.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("sgnRole_PDGF.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("FGF")
pdf("netvisual_agg_FGF.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("netvisual_hmap_FGF.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("sgnRole_FGF.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("CCL")
pdf("netvisual_agg_CCL.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("netvisual_hmap_CCL.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("sgnRole_CCL.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()

pathways.show=c("CXCL")
pdf("netvisual_agg_CXCL.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("netvisual_hmap_CXCL.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("sgnRole_CXCL.pdf")
netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=12,height=2.5,font.size=10)
dev.off()
