##R 
##single cell data of ESCC
##step1: interpretation of the cell cluster
##Jiacheng Dai #daicy0424@gmail.com
##2022/07/05

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

#load data
setwd("01ECatlas/01total_cluster")
matrix_dir = "/home/public/project/singleCell/run/agg77/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir,"barcodes.tsv.gz")
features.path <- paste0(matrix_dir,"features.tsv.gz")
matrix.path <- paste0(matrix_dir,"matrix.mtx.gz")
data.matrix <- readMM(file=matrix.path)
feature.names=read.delim(features.path,header=F,stringsAsFactors=FALSE)
barcode.names=read.delim(barcode.path,header=F,stringsAsFactors=FALSE)
colnames(data.matrix) <- barcode.names$V1
rownames(data.matrix) <- feature.names$V2

#create seurat object
data <- CreateSeuratObject(counts = data.matrix,
                        project = "ESCC",
                        min.cells = 3,
                        min.features = 0)
#filtering
meta = read.csv("escc_metadata.csv",header=T)
meta1 = meta[meta$type!="Blood",]
data = data[,meta1$barcode_id]

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "MT-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "RP[SL]")
pdf("new0418/Vlnplot1.pdf",height=10,width=22)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()
data = subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
# [1]  28900 227671
pdf("new0418/Vlnplot2.pdf",height=10,width=22)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()

#SCTransform, replace the function of NormalizeData(), ScaleData(), FindVariableFeatures()
#data stored at SCT assay
options(future.globals.maxSize = 2000 * 1024^2)
data <- SCTransform(data, vars.to.regress = c("percent.mt","percent.ribo"), verbose = TRUE)
saveRDS(data,"new0418/01ESCCdataset.mt10.rds")

#clustering
data = RunPCA(data)
data = RunUMAP(data, dims = 1:30)
#data = RunTSNE(data, dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(data, resolution = 1) 
saveRDS(data,"new0418/01ESCCdataset.mt10.rds")
meta = read.csv("escc_metadata.csv",header=T,row.names=3)
meta1 = meta[colnames(data),]
write.csv(meta1, "new0418/metadata.mt10.csv")

#Doublet Finder
sample.select = meta1[meta1$Patient_ID == "p055",]
data1 = data[,rownames(sample.select)]
sweep.res.list = paramSweep_v3(data1, PCs = 1:30, sct = TRUE)
sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn = find.pK(sweep.stats)
pK_value = as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))

annotations = data1@meta.data$SCT_snn_res.1
homotypic.prop = modelHomotypic(annotations)
nExp_poi = round(0.075*nrow(data1@meta.data))
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
pN_value = 0.25
pANN_value = paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
data1 <- doubletFinder_v3(data1, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)  
data1 <- doubletFinder_v3(data1, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = TRUE)  
df.classifications = paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
data1@meta.data$Doublet <- data1@meta.data[ ,df.classifications]
head(data1@meta.data) 
table(data1@meta.data$Doublet)
write.csv(data1@meta.data, "new0418/metap055.csv")

fs = list.files()
data.matrix = do.call(rbind, lapply(list.files(), function(x){read.csv(file.path(x), header=T, row.names=1)[,c(8,9,13)]}))
write.csv(data.matrix, "doublet.csv")

doublet = doublet[doublet$Doublet == "Singlet",]
data1 = data[,rownames(doublet)]
saveRDS(data1, "new0418/01ESCCdataset.211750.rds")

#Find cluster biomarkers
cluster.markers = FindAllMarkers(data1,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers,file="01ESCCdataset_Markers.csv")

#plot 
pdf("new0418/01umap.cluster.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()

#
meta = read.csv("new0418/metadata.mt10.csv",header=T,row.names=1)
meta1 = meta[colnames(data),c(4,5,6,10)]
meta1$cluster = data@meta.data$seurat_clusters

meta1.FB = meta1[meta1$cluster %in% c(0,4,5,12,14,23),]
meta1.FB$celltype = "Fibroblast"
meta1.BC = meta1[meta1$cluster %in% c(1,48),]
meta1.BC$celltype = "B cell"
meta1.EC = meta1[meta1$cluster %in% c(6,8,46),]
meta1.EC$celltype = "Endothelial cell"
meta1.epi = meta1[meta1$cluster %in% c(9,11,20,22,37,41),]
meta1.epi$celltype = "Epithelial cell"
meta1.mye = meta1[meta1$cluster %in% c(16,18,19,26,30,44),]
meta1.mye$celltype = "Myeloid cell"
meta1.PC = meta1[meta1$cluster %in% c(24,25,32,33,35,36,39,40,45),]
meta1.PC$celltype = "Plasma cell"
meta1.smc = meta1[meta1$cluster %in% c(13,21,31,43),]
meta1.smc$celltype = "Smooth muscle cell"
meta1.TC = meta1[meta1$cluster %in% c(2,3,7,10,15,28,34),]
meta1.TC$celltype = "T cell"
meta1.2 = meta1[meta1$cluster == 17,]
meta1.2$celltype = "Mast cell"
meta1.3 = meta1[meta1$cluster == 29,]
meta1.3$celltype = "Lymphatic endothelial cell"

meta2 = rbind(meta1.FB, meta1.BC, meta1.EC, meta1.epi, meta1.mye, meta1.PC, meta1.smc, meta1.TC, meta1.2, meta1.3)
meta3.1 = meta2[meta2$type == "Normal",]
meta3.1$tissue = "Distant"
meta3.2 = meta2[meta2$type == "Peri",]
meta3.2$tissue = "Adjacent"
meta3.3 = meta2[meta2$type == "Tumor",]
meta3.3$tissue = "Tumor"
meta3 = rbind(meta3.1, meta3.2, meta3.3) #207115
write.csv(meta3,"new0418/metadata.211.csv")

#plot 
data1 = data[,rownames(meta3)]
data1@meta.data[,c(10:13)] = meta3[,c(4:7)]
pdf("new0418/01umap.cluster1.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "cluster", label=TRUE, pt.size = 0.5)
dev.off()
pdf("new0418/01umap.celltype.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5)
dev.off()
pdf("new0418/01umap.tissuetype.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "tissue", label=TRUE, pt.size = 0.5)
dev.off()
pdf("new0418/01umap.CD45.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "CD45", label=TRUE, pt.size = 0.5)
dev.off()
saveRDS(data1, "new0418/01ESCCdataset.211750.rds")
#
meta = data[[]]
meta1 = meta[meta$tissue=="Tumor",]
meta2 = meta[meta$tissue=="Adjacent",]
meta3 = meta[meta$tissue=="Distant",]
data1 = data[,rownames(meta1)]
data2 = data[,rownames(meta2)]
data3 = data[,rownames(meta3)]
pdf("01umap.dist.pdf",w=8,h=5)
DimPlot(data3, reduction = "umap", group.by = "celltype", label=F, pt.size = 0.5)
dev.off()

#DotPlot
meta$celltype = factor(meta$celltype, levels = c("T cell","Fibroblast","Epithelial cell","Endothelial cell", "Myeloid cell","B cell","Plasma cell","Smooth muscle cell","Mast cell","Lymphatic endothelial cell"))
meta$list = as.numeric(meta$celltype)
for (i in 1:206380) {
     data@active.ident[[i]] = meta$list[i]
}
cd_gene = c("CD3D","CD8A","GZMK","NKG7","CTLA4","DCN","LUM","CFD","MGP","COL1A1","KRT5","KRT6A","KRT13","KRT14","KRT15","PECAM1","CLDN5","CAV1","ACKR1","VWF","C1QA","C1QB","CD163","CXCL8","AIF1","CD74","CD79A","CD37","MS4A1","CXCR4","IGHG1","IGHG2","IGHG3","IGKC","IGHA1","ACTA2","TAGLN","MYL9","MYH11","TPM1","TPSB2","TPSAB1","HPGD","HPGDS","CPA3","CCL21","TFF3","MMRN1","GNG11","CAVIN2")
pdf("DotPlot_all.pdf",h=4,w=20)
DotPlot(data, features=cd_gene) + theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()


#
pdf("new0418/01umap_EC.pdf",w=12,h=5)
FeaturePlot(data,features=c("CDH5","PECAM1"))
dev.off()
pdf("new0418/01umap_SMC.pdf",w=12,h=5)
FeaturePlot(data,features=c("ACTA2","MYH11"))
dev.off()
pdf("new0418/01umap_Mastcell.pdf",w=12,h=5)
FeaturePlot(data,features=c("TPSB2","TPSAB1"))
dev.off()
pdf("new0418/01umap_Fibroblast.pdf",w=12,h=5)
FeaturePlot(data,features=c("LUM","DCN"))
dev.off()
pdf("new0418/01umap_Epi.pdf",w=12,h=10)
FeaturePlot(data,features=c("KRT13","KRT14","KRT6A","EPCAM"))
dev.off()
pdf("new0418/01umap_Tcell.pdf",w=12,h=5)
FeaturePlot(data,features=c("CD8A","GZMK"))
dev.off()
pdf("new0418/01umap_Bcell.pdf",w=12,h=5)
FeaturePlot(data,features=c("CD79A","CD37"))
dev.off()
pdf("new0418/01umap_myeloid.pdf",w=12,h=10)
FeaturePlot(data,features=c("C1QA","C1QB","CXCL8","G0S2"))
dev.off()

pdf("FeaturePlot/01ESCC_Pericyte.umap.pdf",w=6,h=5)
FeaturePlot(data,features="ANPEP")
dev.off()

#cell proportion plot
meta = read.csv("escc_metadata.csv",header=T,row.names=3)
meta = meta[colnames(data),]
meta$seurat_clusters = data$seurat_clusters 
write.csv(meta,"01ESCCdataset_Meta.csv")

library(ggplot2)
library(reshape2)

#1 cluster vs tissue type
rm(list=ls())
options(stringsAsFactors=TRUE)
meta = read.csv("01ESCCdataset_Meta.csv",header=T,row.names=1)
proplist = matrix(NA,nrow=59,ncol=3)
list1 = c(0:58)
j = 1
for (i in list1) {
     meta1 = meta[meta$seurat_clusters==i,]
     proplist[j,] = table(meta1$type)
	   j = j+1
}
colnames(proplist) = c("Distant","Adjacent","Tumor")
rownames(proplist) = paste0("Cluster",c(0:58))
proplist2 = melt(proplist)
colnames(proplist2) = c("seurat_clusters","Type","value")
proplist2$value = as.numeric(proplist2$value)
p <- ggplot(data=proplist2,aes(x=seurat_clusters,y=value,fill=Type))+
     geom_bar(stat="identity")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("propplot/cellabsvalue_cluster_vs_type.pdf",h=5,w=25)
p
dev.off()

p <- ggplot(data=proplist2,aes(x=seurat_clusters,y=value,fill=Type))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("propplot/cellprop_cluster_vs_type.pdf",h=5,w=25)
p
dev.off()

#2 cluster vs patient ID
rm(list=ls())
options(stringsAsFactors=TRUE)
meta = read.csv("01ESCCdataset_Meta.csv",header=T,row.names=1)
proplist = matrix(NA,nrow=59,ncol=20)
list1 = c(0:58)
j = 1
for (i in list1) {
     meta1 = meta[meta$seurat_clusters==i,]
     proplist[j,] = table(meta1$Patient_ID)
	   j = j+1
}
colnames(proplist) = c("p001","p002","p013","p014","p015","p021","p022","p025","p026","p034","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055")
rownames(proplist) = paste0("Cluster",c(0:58))
write.csv(proplist, "propplot/cellmatrix_cluster_vs_patient.csv")

#画简要版的fig大图
rm(list=ls())
options(stringsAsFactors=TRUE)

proplist = data.frame()
cluster = c("T cell", "Mast cell", "B cell", "Epithelial cell", "Myeloid cell", "Endothelial cell", "Lymphatic Endothelial cell", "Fibroblast", "Smooth Muscle Cell", "other")
IDs = c("p001","p002","p013","p014","p015","p021","p022","p025","p026","p034","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055")
type = c("Tumor","Peri","Normal")
j=1
for (a in cluster) {
     for (b in IDs) {
          for (d in type) {
               meta1 = clusters[which(clusters$celltype==a & clusters$Patient_ID==b & clusters$type==d),]
               proplist[j,1:4] = c(a,b,d,nrow(meta1))
               j=j+1
          }
     }
}
colnames(proplist) = c("CellType","PatientID","TissueType","CellNumbers")
write.csv(proplist, "propplot/01ESCCdataset_cellprop_Cluster_ID_Type.csv")

#按分成的10类来画Dimplot图
data = data[,rownames(clusters)]
data@meta.data$CellType = clusters$celltype
pdf("01ESCC_umap_celltype.pdf",w=12,h=10)
DimPlot(data, reduction="umap", label=TRUE, label.size=4, group.by="CellType", pt.size=.1)
dev.off()
data@meta.data$TissueType = clusters$type
pdf("01ESCC_umap_tissuetype.pdf",w=12,h=10)
DimPlot(data, reduction="umap", group.by="TissueType", pt.size=.1)
dev.off()

#按10类来画heatmap图
data = readRDS("new0418/01ESCCdataset.211750.rds")
meta = read.csv("new0418/metadata.211.csv",header=T,row.names=1)
cluster1 = meta[meta$celltype == "Fibroblast",]
cluster1$cluster = 1
cluster2 = meta[meta$celltype == "B cell",]
cluster2$cluster = 2
cluster3 = meta[meta$celltype == "Endothelial cell",]
cluster3$cluster = 3 
cluster4 = meta[meta$celltype == "Epithelial cell",]
cluster4$cluster = 4
cluster5 = meta[meta$celltype == "Myeloid cell",]
cluster5$cluster = 5
cluster6 = meta[meta$celltype == "Plasma cell",]
cluster6$cluster = 6
cluster7 = meta[meta$celltype == "Smooth muscle cell",]
cluster7$cluster = 7
cluster8 = meta[meta$celltype == "T cell",]
cluster8$cluster = 8
cluster9 = meta[meta$celltype == "Mast cell",]
cluster9$cluster = 9
cluster10 = meta[meta$celltype == "Lymphatic endothelial cell",]
cluster10$cluster = 10
cluster11 = meta[meta$celltype == "Fibroblast.independent",]
cluster11$cluster = 11

cluster = rbind(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8, cluster9, cluster10)
data = data[,rownames(cluster)]
for (i in 1:206380) {
     data@active.ident[[i]] = cluster[i,10]
}
data@active.ident = factor(data@active.ident, levels = c(1,2,3,4,5,6,7,8,9,10))
levels(data@active.ident)
data.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(data.markers, "01ESCCdataset_marker_10cc.csv")
top = data.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
set.seed(424)
subdata = subset(data, downsample = 400)
plot = DoHeatmap(subdata, features=top$gene) + NoLegend()
pdf("new0418/Heatmap.pdf",h=10,w=12)
plot
dev.off()


proplist$CellNumbers = as.numeric(proplist$CellNumbers)
proplist1 = proplist[proplist$TissueType=="Normal",]
p <- ggplot(data=proplist1,aes(x=PatientID,y=CellNumbers,fill=CellType))+
     geom_bar(stat="identity")+
     theme_minimal()+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("propplot/cellabsvalue_Patient_distant.pdf",h=5,w=10)
p
dev.off()

#饼图
data = c(54278, 42386, 51068, 95866, 4995, 15522, 38236, 5781, 37743, 92138)
names = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
cols = c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3","#FF62BC")
pdf("propplot/pieplot_all.pdf",h=5,w=8)
pie(data, labels=names, col=cols)
dev.off()

data = c(19758, 11252, 32282, 23946, 600, 3009, 16921, 735, 9493, 42028)
names = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
cols = c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3","#FF62BC")
pdf("propplot/pieplot_tumor.pdf",h=5,w=8)
pie(data, labels=names, col=cols)
dev.off()

data = c(18378, 12885, 6301, 22123, 1406, 4364, 8196, 1040, 9793, 18757)
names = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
cols = c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3","#FF62BC")
pdf("propplot/pieplot_adjacent.pdf",h=5,w=8)
pie(data, labels=names, col=cols)
dev.off()

data = c(16142, 17149, 12485, 49797, 2989, 8149, 13119, 4006, 18457, 31353)
names = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
cols = c("#F8766D","#D89000","#A3A500","#39B600","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3","#FF62BC")
pdf("propplot/pieplot_distant.pdf",h=5,w=8)
pie(data, labels=names, col=cols)
dev.off()

#按10类来画cellprop和cellabsvalue
proplist = matrix(NA,nrow=10,ncol=3)
list1 = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
j = 1
for (i in list1) {
     meta1 = clusters[clusters$celltype==i,]
     proplist[j,] = table(meta1$type)
	   j = j+1
}
colnames(proplist) = c("Distant","Adjacent","Tumor")
rownames(proplist) = c("B cell","Endothelial cell","Epithelial cell","Fibroblast","Lymphatic Endothelial cell","Mast cell","Myeloid cell","other","Smooth Muscle Cell","T cell")
proplist2 = melt(proplist)
colnames(proplist2) = c("CellType","TissueType","value")
proplist2$value = as.numeric(proplist2$value)
p1 <- ggplot(data=proplist2,aes(x=CellType,y=value,fill=TissueType))+
     geom_bar(stat="identity")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("propplot/cellabsvalue_cluster_vs_type_10cc.pdf",h=5,w=8)
p1
dev.off()

p2 <- ggplot(data=proplist2,aes(x=CellType,y=value,fill=TissueType))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("propplot/cellprop_cluster_vs_type_10cc.pdf",h=5,w=8)
p2
dev.off()

p2 <- ggplot(data=proplist,aes(x=PatientID,y=CellNumbers,fill=CellType))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#E6194B","#3CB44B","#FFE119","#4363D8","#F58231","#911EB4","#42D4F4","#BFEF45","#FABED4","#F032E6"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#average expression
data = readRDS("01ESCCdataset.211750.rds")
meta = data[[]]
meta = meta[meta$celltype != "Fibroblast.independent",]
data = data[,rownames(meta)]
metadata=read.csv("metadata.211.csv",header=T,row.names=1)
metadata=metadata[rownames(data@meta.data),]
data@meta.data$Patient_ID = metadata$Patient_ID
meta$celltype = factor(meta$celltype, levels = c("T cell","Fibroblast","Epithelial cell","Endothelial cell", "Myeloid cell","B cell","Plasma cell","Smooth muscle cell","Mast cell","Lymphatic endothelial cell"))
meta$list = as.numeric(meta$celltype)
for (i in 1:206380) {
     data@active.ident[[i]] = meta$list[i]
}
data@active.ident = factor(data@active.ident, levels = c(1,2,3,4,5,6,7,8,9,10))

#
pid = c("p026","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055")
tis = c("Tumor","Adjacent","Distant")
meta2 = data[[]]

meta3 = meta2[meta2$Patient_ID=="p001" & meta2$tissue=="Tumor",]
data1 = data[,rownames(meta3)]
test = AverageExpression(data1)

for (i in 1:11) {
     for (j in 1:3) {
          meta3 = meta2[meta2$Patient_ID==pid[i] & meta2$tissue==tis[j],]
          data1 = data[,rownames(meta3)]
          test = AverageExpression(data1)
          saveRDS(test, file=paste0("hist/AvExp.",pid[i],tis[j],".rds"))
     }
}
