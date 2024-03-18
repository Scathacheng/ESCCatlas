#R
#cellchat

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(CellChat) #v1.0.0
library(ggalluvial)
library(ggplot2)

#1 combine data

#epithelial
meta = read.csv("02epithelial/metadata.res1.3.csv",header=T,row.names=1)
meta1 = meta[,c("patient","tissue","major","epi_celltype")]
colnames(meta1) = c("patient","tissue","celltype","subcelltype")

#T cell
meta = read.csv("03tcell/metadata.res1.4.csv",header=T,row.names=1)
meta2 = meta[,c("patient","tissue","major","tcell_celltype")]
colnames(meta2) = c("patient","tissue","celltype","subcelltype")

#B cell
meta = read.csv("04bcell/metadata.res1.1.csv",header=T,row.names=1)
meta3 = meta[,c("patient","tissue","major","bcell_celltype")]
colnames(meta3) = c("patient","tissue","celltype","subcelltype")

#Fibroblast
meta = read.csv("05fibroblast/metadata.res1.3.csv",header=T,row.names=1)
meta4 = meta[,c("patient","tissue","major","fb_celltype")]
colnames(meta4) = c("patient","tissue","celltype","subcelltype")

#myeloid cell
meta = read.csv("06myeloidcell/metadata.res1.3.csv",header=T,row.names=1)
meta5 = meta[,c("patient","tissue","major","myeloid_subtype")]
colnames(meta5) = c("patient","tissue","celltype","subcelltype")

#endothelial cell
meta = read.csv("total.16580.csv",header=T,row.names=1)
meta6 = meta[,c("Patient_ID","tissue","orig.ident","celltype")]
colnames(meta6) = c("patient","tissue","celltype","subcelltype")

meta = rbind(meta1, meta2, meta3, meta4, meta5, meta6)
#185714
write.csv(meta, "08cci/metadata.185714.csv")



#2 load data
rm(list=ls())
meta = read.csv("metadata.118364.csv",header=T,row.names=1)
data = readRDS("../new0418/01ESCCdataset.211750.rds")
data1 = data[,rownames(meta)]
matrix1 = data1@assays$SCT@data

meta1 = meta[meta$tissue =="Tumor",]
data.input = matrix1[,rownames(meta1)]
cellchat <- createCellChat(object = data.input, meta = meta1, group.by = "celltype")

#add cell information
cellchat <- addMeta(cellchat, meta = meta1)
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

df.net <- subsetCommunication(cellchat, sources.use = c(1:7), targets.use = c(1:7))
#infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#calculate the network
cellchat <- aggregateNet(cellchat,sources.use = c(1:7), targets.use = c(1:7), remove.isolate=FALSE)

### visualization of cell-cell communication network
pdf("bubble_Tumor.pdf",h=20,w=15)
netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(1:7), remove.isolate = FALSE)
dev.off()
pdf("chord_Tumor.pdf",h=10,w=10)
netVisual_chord_gene(cellchat, sources.use = c(1:7), targets.use = c(1:7), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

#cellphoneDB
rm(list=ls())
meta = read.csv("08cci/metadata.185714.csv",header=T,row.names=1)
data = readRDS("new1228/01Alldata.harmony.rds")
overlap = intersect(rownames(meta),colnames(data))
meta = meta[overlap,]
data1 = data[,overlap]
matrix1 = data1@assays$SCT@data

meta1 = meta[meta$tissue =="Tumor",]
data.input = matrix1[,rownames(meta1)]
gene = rownames(data.input)
rt1 = cbind(gene, as.matrix(data.input))
write.table(rt1, "08cci/cellphonedb_count_Tumor.txt", sep = "\t", quote = F, row.names = F, col.names = T)
test = meta1[,c("patient","subcelltype")]
colnames(test)=c("sampleID","celltype")
test[,1] = rownames(meta1)
write.table(test, "08cci/cellphonedb_meta_Tumor.txt", sep = "\t", quote = F, row.names = F)

##4 cellphoneDB
source anaconda2/bin/activate
source activate cellphonedb3
cellphonedb method statistical_analysis cellphonedb_meta_Dist.txt cellphonedb_count_Dist.txt --threads 32 --counts-data gene_name
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot cellphonedb_meta_Tumor.txt