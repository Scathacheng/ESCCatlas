##R 
##single cell data of ESCC
##step2: EC sub-cluster
##Jiacheng Dai #daicy0424@gmail.com
##2022/07/10

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(DoubletFinder)

#load data 
setwd("01ECatlas/02EC_cluster")
data = readRDS("../01total_cluster/new0418/01ESCCdataset.211750.rds")
meta = read.csv("01total_cluster/new0418/metadata.211.csv",header=T,row.names=1)
meta.ec = meta[meta$celltype == "Endothelial cell",] #18537
data.matrix = data@assays$RNA@counts[,rownames(meta.ec)]
#dim(data.matrix) 28900 46281

data1 <- CreateSeuratObject(counts = data.matrix,
                        project = "endothelial cells",
                        min.cells = 3,
                        min.features = 0)
options(future.globals.maxSize = 2000 * 1024^2)
data1 <- SCTransform(data1, verbose = T)

#clustering
data1 = RunPCA(data1)
data1 = RunUMAP(data1, dims = 1:30)
#data.ec = RunTSNE(data.ec, dims = 1:30)
data1 = FindNeighbors(data1, dims = 1:30) 
data1 = FindClusters(data1, resolution = 1) 
saveRDS(data1,"new0423/02ECdataset.sct.rds")

#Find cluster biomarkers
cluster.markers = FindAllMarkers(data1,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers,file="new0423/EC_Markers.csv")

#plot 
pdf("new0423/02umap.SCT1.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "SCT_snn_res.1", label=TRUE, pt.size = 0.5)
dev.off()

pdf("new0423/02umap.SCT05.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "SCT_snn_res.0.5", label=TRUE, pt.size = 0.5)
dev.off()

#
data = readRDS("02EC_cluster/new0423/02ECdataset.sct.rds")
meta = read.csv("01total_cluster/new0418/metadata.211.csv", header=T, row.names=1)
meta = meta[meta$celltype == "Endothelial cell",]
meta$seurat_clusters = data@meta.data$seurat_clusters
meta0 = meta[meta$seurat_clusters == 0,]; meta0$celltype = "EC_IGKC+"
meta1 = meta[meta$seurat_clusters == 1,]; meta1$celltype = "Artery_DEPP1+"
meta2 = meta[meta$seurat_clusters == 2,]; meta2$celltype = "Vein_HLA+"
meta3 = meta[meta$seurat_clusters == 3,]; meta3$celltype = "stalk cell_RGCC+"
meta4 = meta[meta$seurat_clusters == 4,]; meta4$celltype = "Capillary_FABP5+"
meta5 = meta[meta$seurat_clusters == 5,]; meta5$celltype = "Vein_SELE+"
meta6 = meta[meta$seurat_clusters == 6,]; meta6$celltype = "tip cell_COL4A1+"
meta7 = meta[meta$seurat_clusters == 7,]; meta7$celltype = "Vein_IL6+"
meta8 = meta[meta$seurat_clusters == 8,]; meta8$celltype = "Capillary_SOCS3+"
meta9 = meta[meta$seurat_clusters == 9,]; meta9$celltype = "Vein_CPE+"
meta10 = meta[meta$seurat_clusters == 10,]; meta10$celltype = "Vein_ISG15+"
meta11 = meta[meta$seurat_clusters == 11,]; meta11$celltype = "PCV_ACTG1+"
meta12 = meta[meta$seurat_clusters == 12,]; meta12$celltype = "Capillary_CA4+"
meta13 = meta[meta$seurat_clusters == 13,]; meta13$celltype = "Artery_GAS6+"
meta14 = meta[meta$seurat_clusters == 14,]; meta14$celltype = "EC_HSP+"
meta18 = meta[meta$seurat_clusters == 18,]; meta18$celltype = "Artery_UNC5B+"
total = rbind(meta1, meta13, meta18, meta0, meta2, meta5, meta7, meta9, meta10, meta11, meta4, meta8, meta12, meta6, meta3, meta14) #16580
data1 = data[,rownames(total)]
data1@meta.data$celltype = total$celltype
data1@meta.data$tissue = total$tissue
saveRDS(data1, "02EC_cluster/new0423/02ECdataset.16580.rds")

#plot
pdf("new0423/02umap.cluster1.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()
pdf("new0423/02umap.celltype.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("new0423/02umap.tissuetype.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "tissue", label=TRUE, pt.size = 0.5)
dev.off()

#FeaturePlot
pdf("new0423/02umap.featureHSPG2.pdf",w=6,h=5)
FeaturePlot(data1,features=c("HSPG2"))
dev.off()
pdf("new0423/02umap.featureINSR.pdf",w=6,h=5)
FeaturePlot(data1,features=c("INSR"))
dev.off()
pdf("new0423/02umap.featureCOL4A1.pdf",w=6,h=5)
FeaturePlot(data1,features=c("COL4A1"))
dev.off()
pdf("new0423/02umap.featureGJA5.pdf",w=6,h=5)
FeaturePlot(data1,features=c("GJA5"))
dev.off()
pdf("new0423/02umap.featureACKR1.pdf",w=6,h=5)
FeaturePlot(data1,features=c("ACKR1"))
dev.off()
pdf("new0423/02umap.featureCA4.pdf",w=6,h=5)
FeaturePlot(data1,features=c("CA4"))
dev.off()


#Dotplot
library(Seurat)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggplot2)
library(patchwork)
#rm(list=ls())
#artery
cd_gene<-c("SOX17","EFNB2","GJA5","JAG1","SEMA3G","DEPP1","GAS6","FBLN2","LTBP4","FN1","SERPINE2","UNC5B","CXCL12","DLL4","NOTCH4","CTNNB1","ENG")
#vein
cd_gene<-c("ACKR1","VWF","HLA-DQA1","HLA-DRB5","SELE","SELP","PRCP","IL6","NNMT","ACKR3","CPE","IL1R1","AKR1C1","ISG15","IFI44L","LY6E","ACTG1","ANGPT2","CAV1")
#capillary
cd_gene<-c("CA4","CAV1","CAV2","FABP4","FABP5","GSN","SOCS3","CXCL12","SPRY1","CXCL2","CAVIN2","ANXA3")
#TEC
cd_gene<-c("COL4A1","COL4A2","SERPINH1","MMP2","SPARC","INSR","HSPG2","NOTCH4","PLVAP","VWA1","IGFBP5","APLNR","PMEPA1","KDR","FLT1","TIE1")
#unknown
cd_gene<-c("IGKC","HSPA1A","HSPA1B","HSP90AA1","HSPD1","HSPE1","CCL2","ATF3","CYR61","UBC")

#plot
data.usage = DotPlot(data1,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(1,13,18,2,5,7,9,10,11,4,8,12,6,3,0,14))
levels(data.usage$id) = c("Artery_DEPP1","Artery_GAS6","Artery_UNC5B","Vein_HLA","Vein_SELE","Vein_IL6","Vein_CPE","Vein_ISG15","PCV_ACTG1","Capillary_FABP5","Capillary_SOCS3","Capillary_CA4","tip cell_COL4A1","stalk EC_RGCC","EC_IGKC","EC_HSP")

pdf("new0423/dotplot_TEC2.pdf",h=6,w=8)
data.usage %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5)) +
theme(axis.title.x = element_blank()) +
ylab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(0,3))
dev.off()


##cell proportion plot
meta = read.csv("../01total_cluster/01ESCCdataset_Meta.csv",header=T,row.names=1)
meta = meta[colnames(data),]
meta$seurat_clusters_EC = data$seurat_clusters 
meta$Doublet = data$Doublet
write.csv(meta,"02ECdataset_Meta.csv")

library(ggplot2)
library(reshape2)

#1 cluster vs tissue type
rm(list=ls())
options(stringsAsFactors=TRUE)
meta = read.csv("02ECdataset_Meta.csv",header=T,row.names=1)
proplist = matrix(NA,nrow=25,ncol=3)
list1 = c(0:24)
j = 1
for (i in list1) {
     meta1 = meta[meta$seurat_clusters_EC==i,]
     proplist[j,] = table(meta1$type)
	   j = j+1
}
colnames(proplist) = c("Distant","Adjacent","Tumor")
rownames(proplist) = paste0("ECcluster",c(0:24))
proplist2 = melt(proplist)
colnames(proplist2) = c("seurat_clusters","Type","value")
proplist2$value = as.numeric(proplist2$value)
p <- ggplot(data=proplist2,aes(x=seurat_clusters,y=value,fill=Type))+
     geom_bar(stat="identity")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("ECpropplot/cellabsvalue_cluster_vs_type.pdf",h=5,w=7)
p
dev.off()

p <- ggplot(data=proplist2,aes(x=seurat_clusters,y=value,fill=Type))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#008000","#0000FF","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("ECpropplot/cellprop_cluster_vs_type.pdf",h=5,w=7)
p
dev.off()

#2 cluster vs patient ID
rm(list=ls())
options(stringsAsFactors=TRUE)
meta = read.csv("02ECdataset_Meta.csv",header=T,row.names=1)
meta = meta[meta$Doublet=="Singlet",]

proplist = data.frame()
clusters = c(0:24)
IDs = c("p001","p002","p013","p014","p015","p021","p022","p025","p026","p034","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055")
type = c("Tumor","Peri","Normal")
j=1
for (a in clusters) {
     for (b in IDs) {
          for (d in type) {
               meta1 = meta[which(meta$seurat_clusters_EC==a & meta$Patient_ID==b & meta$type==d),]
               proplist[j,1:4] = c(a,b,d,nrow(meta1))
               j=j+1
          }
     }
}
colnames(proplist) = c("ECcluster","PatientID","TissueType","CellNumbers")
write.csv(proplist,"ECpropplot/02ECdataset_cellprop_Cluster_ID_Type.csv")

#statistical analysis
#1 normalize cell number
proplist$CellNumbers = as.numeric(proplist$CellNumbers)
clusters = c(0:24)
IDs = c("p001","p002","p013","p014","p015","p021","p022","p025","p026","p034","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055")
type = c("Tumor","Peri","Normal")
for (a in type) {
     for (b in IDs) {
          proplist1 = proplist[which(proplist$PatientID==b & proplist$TissueType==a),]
          i = sum(proplist1$CellNumbers)
          for (d in clusters) {
               proplist[proplist$PatientID==b & proplist$TissueType==a & proplist$ECcluster==d,5] = proplist[proplist$PatientID==b & proplist$TissueType==a & proplist$ECcluster==d,4]/i
          }
     }
}
colnames(proplist)[5] = "normCellProp"

#2 t test and wilcox test
statistic = data.frame()
j=1
for (a in clusters) {
     proplist_C = proplist[proplist$ECcluster==a,]
     proplist_C_T = proplist_C[proplist_C$TissueType=="Tumor",]
     proplist_C_D = proplist_C[proplist_C$TissueType=="Normal",]
     d = wilcox.test(proplist_C_T$normCellProp,proplist_C_D$normCellProp)
     e = t.test(proplist_C_T$normCellProp,proplist_C_D$normCellProp,paired=TRUE, alternative="two.sided")
     statistic[j,1:6] = c(a, "TumorvsDistant", d$statistic, d$p.value, e$statistic, e$p.value)
     j=j+1
}
colnames(statistic) = c("ECcluster","Compare","Wilcox.stc", "Wilcox.p", "Ttest.stc", "Ttest.p")
statistic$Wilcox.p = as.numeric(statistic$Wilcox.p)
statistic$Wilcox.fdr = p.adjust(statistic$Wilcox.p,method="fdr")
statistic$Ttest.p = as.numeric(statistic$Ttest.p)
statistic$Ttest.fdr = p.adjust(statistic$Ttest.p,method="fdr")
write.csv(statistic, "ECpropplot/02ECdataset_cellprop_statistic.csv")

#3 histogram
proplist$CellNumbers = as.numeric(proplist$CellNumbers)
proplist$ECcluster = as.factor(proplist$ECcluster)
levels(proplist$ECcluster) = c("HLA+ vein (C0)","Qsc EC (C1)","CPE+ vein (C10)","ANGPT2+ LEC (C11)","ISG15+ vein (C12)","MGP+ artery (C13)","HSP+ EC (C14)","KRT+ EC (C15)","ACTG1+ vein (C16)","Fibroblast-like (C17)","Immune cell-like (C18)","APLNR+ TEC (C19)","RGCC+ EC (C2)","SMC-like (C20)","B cell-like (C21)","SELE+ vein (C22)","Macrophage-like (C23)","TOP2A+ TEC (C24)","COL4A1+ TEC (C3)","CCL2+ vein (C4)","JAG1+ artery (C5)","CAV1+ capillary (C6)","LYVE1+ LEC (C7)","FABP4+ EC (C8)","SOCS3+ EC (C9)")
proplist$TissueType = as.factor(proplist$TissueType)
levels(proplist$TissueType) = c("Distant","Adjacent","Tumor")

celltypelist = c("HLA+ vein (C0)","Qsc EC (C1)","CPE+ vein (C10)","ANGPT2+ LEC (C11)","ISG15+ vein (C12)","MGP+ artery (C13)","HSP+ EC (C14)","KRT+ EC (C15)","ACTG1+ vein (C16)","Fibroblast-like (C17)","Immune cell-like (C18)","APLNR+ TEC (C19)","RGCC+ EC (C2)","SMC-like (C20)","B cell-like (C21)","SELE+ vein (C22)","Macrophage-like (C23)","TOP2A+ TEC (C24)","COL4A1+ TEC (C3)","CCL2+ vein (C4)","JAG1+ artery (C5)","CAV1+ capillary (C6)","LYVE1+ LEC (C7)","FABP4+ EC (C8)","SOCS3+ EC (C9)")

for (i in 1:25) {
     proplist1 = proplist[proplist$ECcluster==celltypelist[i],]
     p <- ggplot(data=proplist1,aes(x=PatientID,y=CellNumbers,fill=TissueType))+
     geom_bar(stat="identity")+
     theme_minimal()+
     scale_fill_manual(values=c("#0000FF","#008000","#ff0000"))+
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
     pdf(paste0("cellprop_",celltypelist[i],".pdf"),h=5,w=8)
     p
     dev.off()
}


#histogram for signature genes' pathways
library(clusterProfiler)
library(org.Hs.eg.db)
data=read.csv("02ECdataset_Markers.csv",header=T,row.names=1)


#TEC C3 and C19
####################################################################################################################
c3=data[data$cluster==3,]
go3=enrichGO(c3$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go3$qvalue[1:5])
pdf("GO_C3.pdf",h=4,w=5)
barplot(value, col="#CF9400", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go3, "GO_C3.csv")

c19=data[data$cluster==19,]
go19=enrichGO(c19$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go19$qvalue[1:5])
pdf("GO_C19.pdf",h=4,w=5)
barplot(value, col="#CF78FF", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go19, "GO_C19.csv")

c22=data[data$cluster==22,]
go22=enrichGO(c22$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go22$qvalue[c(1,4,5,11,21)])
pdf("GO_C22_2.pdf",h=4,w=5)
barplot(value, col="#FF61C9", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go22, "GO_C22.csv")

c24=data[data$cluster==24,]
go24=enrichGO(c24$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go24$qvalue[1:5])
pdf("GO_C24.pdf",h=4,w=5)
barplot(value, col="#FF6C91", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go24, "GO_C24.csv")
####################################################################################################################


#vein C0 C4 C10 C12 C14 C16
####################################################################################################################
c0=data[data$cluster==0,]
go0=enrichGO(c0$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go0$qvalue[c(2,8,11,12,15)])
pdf("GO_C0.pdf",h=4,w=5)
barplot(value, col="#F8766D", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go0, "GO_C0.csv")

c4=data[data$cluster==4,]
go4=enrichGO(c4$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go4$qvalue[c(2,3,10,11,12)])
pdf("GO_C4_2.pdf",h=4,w=5)
barplot(value, col="#BB9D00", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go4, "GO_C4.csv")

c10=data[data$cluster==10,]
go10=enrichGO(c10$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go10$qvalue[c(1,4,5,6,13)])
pdf("GO_C10_2.pdf",h=4,w=5)
barplot(value, col="#00BF7D", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go10, "GO_C10.csv")

c12=data[data$cluster==12,]
go12=enrichGO(c12$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go12$qvalue[c(3,13,26,28,29)])
pdf("GO_C12_2.pdf",h=4,w=5)
barplot(value, col="#00C0B8", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go12, "GO_C12.csv")

c14=data[data$cluster==14,]
go14=enrichGO(c14$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go14$qvalue[1:5])
pdf("GO_C14.pdf",h=4,w=5)
barplot(value, col="#00B8E5", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go14, "GO_C14.csv")

c16=data[data$cluster==16,]
go16=enrichGO(c16$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go16$qvalue[1:5])
pdf("GO_C16.pdf",h=4,w=5)
barplot(value, col="#03A5FF", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go16, "GO_C16.csv")
####################################################################################################################

#ProfEC, C2, C8, C9
####################################################################################################################
c2=data[data$cluster==2,]
go2=enrichGO(c2$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go2$qvalue[1:5])
pdf("GO_C2.pdf",h=4,w=5)
barplot(value, col="#E08B00", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go2, "GO_C2.csv")

c8=data[data$cluster==8,]
go8=enrichGO(c8$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go8$qvalue[1:5])
pdf("GO_C8.pdf",h=4,w=5)
barplot(value, col="#00B81F", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go8, "GO_C8.csv")

c9=data[data$cluster==9,]
go9=enrichGO(c9$gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
value = -log10(go9$qvalue[1:5])
pdf("GO_C9.pdf",h=4,w=5)
barplot(value, col="#00BC59", horiz=TRUE, cex.names=0.8)
dev.off()
write.csv(go9, "GO_C9.csv")
####################################################################################################################

#clustering robustness evaluation
library(Seurat) # version 4.3.0
data = readRDS("02ECdataset.sct.rds")

#0
dir.create("raw")
write.csv(data@meta.data,"raw/clusters.csv")
pdf("raw/UMAP_raw.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, label.size=5)
dev.off()

#1 
data = FindNeighbors(data, annoy.metric="euclidean", dims=1:30) #28
data = FindClusters(data, resolution=1)
dir.create("euclidean")
write.csv(data@meta.data,"euclidean/clusters.csv")
pdf("euclidean/UMAP_euclidean.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, label.size=5)
dev.off()

#2
data = FindNeighbors(data, annoy.metric="cosine", dims=1:30) #30
data = FindClusters(data, resolution=1)
dir.create("cosine")
write.csv(data@meta.data,"cosine/clusters.csv")
pdf("cosine/UMAP_euclidean.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, label.size=5)
dev.off()

#3
data = FindNeighbors(data, annoy.metric="manhattan", dims=1:30) #28
data = FindClusters(data, resolution=1)
dir.create("manhattan")
write.csv(data@meta.data,"manhattan/clusters.csv")
pdf("manhattan/UMAP_manhattan.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, label.size=5)
dev.off()

#4
data = FindNeighbors(data, annoy.metric="hamming", dims=1:30) #27
data = FindClusters(data, resolution=1)
dir.create("hamming")
write.csv(data@meta.data,"hamming/clusters.csv")
pdf("hamming/UMAP_hamming.pdf",w=12,h=10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, label.size=5)
dev.off()

#combine the clusters.csv file

#plot 1: default parameter and euclidean
library(WGCNA)
rm(list=ls())
clusters = read.csv("clusters.csv",header=T,row.names=1)
clusters.raw = c("JAG1+ artery","MGP+ artery","HLA+ vein","CCL2+ vein","CPE+ vein","ISG15+ vein","ACTG1+ vein","SELE+ vein","CAV1+ capillary","COL4A1+ TEC","APLNR+ TEC","TOP2A+ TEC","LYVE1+ LEC","ANGPT2+ LEC","Qsc EC","RGCC+ EC","FABP4+ EC","SOCS3+ EC","HSP+ EC","KRT+ EC","Fibroblast-like","Immune cell-like","SMC-like","B cell-like","Macrophage-like")
clusters.euc = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27")
n.raw = length(clusters.raw)
n.euc = length(clusters.euc)
pTable = matrix(0, nrow = n.raw, ncol = n.euc)
countbl = matrix(0, nrow = n.raw, ncol = n.euc)
for (fmod in 1:n.raw) 
  for (cmod in 1:n.euc)
  {
    femMembers = (clusters$raw_para == clusters.raw[fmod])
    consMembers = (clusters$euclidean == clusters.euc[cmod])
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    countbl[fmod, cmod] = sum(clusters$raw_para == clusters.raw[fmod] & clusters$euclidean == clusters.euc[cmod])
  }
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50

totals.raw = apply(countbl,1,sum)
totals.euc = apply(countbl,2,sum)
pdf(file="euclidean/consersus.pdf",h=8,w=14)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)

labeledHeatmap(Matrix = pTable,
yLabels = paste(" ", clusters.raw),
xLabels = paste(" ", clusters.euc),
colorLabels = TRUE,
ySymbols = paste("ESCC_", clusters.raw, ": ", totals.raw, sep=""),
xSymbols = paste("EAC_", clusters.euc, ": ", totals.euc, sep=""),
textMatrix = countbl,
colors = greenWhiteRed(100)[50:100],
main = "Cell cluster overlap of default parameter and euclidean",
cex.text = 1.0,
cex.lab = 1.0,
setStdMargins = FALSE)
dev.off()


#plot 2: default parameter and cosine
rm(list=ls())
clusters = read.csv("clusters.csv",header=T,row.names=1)
clusters.raw = c("JAG1+ artery","MGP+ artery","HLA+ vein","CCL2+ vein","CPE+ vein","ISG15+ vein","ACTG1+ vein","SELE+ vein","CAV1+ capillary","COL4A1+ TEC","APLNR+ TEC","TOP2A+ TEC","LYVE1+ LEC","ANGPT2+ LEC","Qsc EC","RGCC+ EC","FABP4+ EC","SOCS3+ EC","HSP+ EC","KRT+ EC","Fibroblast-like","Immune cell-like","SMC-like","B cell-like","Macrophage-like")
clusters.cosine = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28")
n.raw = length(clusters.raw)
n.cosine = length(clusters.cosine)
pTable = matrix(0, nrow = n.raw, ncol = n.cosine)
countbl = matrix(0, nrow = n.raw, ncol = n.cosine)
for (fmod in 1:n.raw) 
  for (cmod in 1:n.cosine)
  {
    femMembers = (clusters$raw_para == clusters.raw[fmod])
    consMembers = (clusters$cosine == clusters.cosine[cmod])
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    countbl[fmod, cmod] = sum(clusters$raw_para == clusters.raw[fmod] & clusters$cosine == clusters.cosine[cmod])
  }
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50

totals.raw = apply(countbl,1,sum)
totals.cosine = apply(countbl,2,sum)
pdf(file="cosine/consersus.pdf",h=8,w=14)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)

labeledHeatmap(Matrix = pTable,
yLabels = paste(" ", clusters.raw),
xLabels = paste(" ", clusters.cosine),
colorLabels = TRUE,
ySymbols = paste("ESCC_", clusters.raw, ": ", totals.raw, sep=""),
xSymbols = paste("EAC_", clusters.cosine, ": ", totals.cosine, sep=""),
textMatrix = countbl,
colors = greenWhiteRed(100)[50:100],
main = "Cell cluster overlap of default parameter and cosine",
cex.text = 1.0,
cex.lab = 1.0,
setStdMargins = FALSE)
dev.off()


#plot 3: default parameter and manhattan
rm(list=ls())
clusters = read.csv("clusters.csv",header=T,row.names=1)
clusters.raw = c("JAG1+ artery","MGP+ artery","HLA+ vein","CCL2+ vein","CPE+ vein","ISG15+ vein","ACTG1+ vein","SELE+ vein","CAV1+ capillary","COL4A1+ TEC","APLNR+ TEC","TOP2A+ TEC","LYVE1+ LEC","ANGPT2+ LEC","Qsc EC","RGCC+ EC","FABP4+ EC","SOCS3+ EC","HSP+ EC","KRT+ EC","Fibroblast-like","Immune cell-like","SMC-like","B cell-like","Macrophage-like")
clusters.manhattan = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27")
n.raw = length(clusters.raw)
n.manhattan = length(clusters.manhattan)
pTable = matrix(0, nrow = n.raw, ncol = n.manhattan)
countbl = matrix(0, nrow = n.raw, ncol = n.manhattan)
for (fmod in 1:n.raw) 
  for (cmod in 1:n.manhattan)
  {
    femMembers = (clusters$raw_para == clusters.raw[fmod])
    consMembers = (clusters$manhattan == clusters.manhattan[cmod])
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    countbl[fmod, cmod] = sum(clusters$raw_para == clusters.raw[fmod] & clusters$manhattan == clusters.manhattan[cmod])
  }
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50

totals.raw = apply(countbl,1,sum)
totals.manhattan = apply(countbl,2,sum)
pdf(file="manhattan/consersus.pdf",h=8,w=14)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)

labeledHeatmap(Matrix = pTable,
yLabels = paste(" ", clusters.raw),
xLabels = paste(" ", clusters.manhattan),
colorLabels = TRUE,
ySymbols = paste("ESCC_", clusters.raw, ": ", totals.raw, sep=""),
xSymbols = paste("EAC_", clusters.manhattan, ": ", totals.manhattan, sep=""),
textMatrix = countbl,
colors = greenWhiteRed(100)[50:100],
main = "Cell cluster overlap of default parameter and manhattan",
cex.text = 1.0,
cex.lab = 1.0,
setStdMargins = FALSE)
dev.off()


#plot 4: default parameter and hamming
rm(list=ls())
clusters = read.csv("clusters.csv",header=T,row.names=1)
clusters.raw = c("JAG1+ artery","MGP+ artery","HLA+ vein","CCL2+ vein","CPE+ vein","ISG15+ vein","ACTG1+ vein","SELE+ vein","CAV1+ capillary","COL4A1+ TEC","APLNR+ TEC","TOP2A+ TEC","LYVE1+ LEC","ANGPT2+ LEC","Qsc EC","RGCC+ EC","FABP4+ EC","SOCS3+ EC","HSP+ EC","KRT+ EC","Fibroblast-like","Immune cell-like","SMC-like","B cell-like","Macrophage-like")
clusters.hamming = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26")
n.raw = length(clusters.raw)
n.hamming = length(clusters.hamming)
pTable = matrix(0, nrow = n.raw, ncol = n.hamming)
countbl = matrix(0, nrow = n.raw, ncol = n.hamming)
for (fmod in 1:n.raw) 
  for (cmod in 1:n.hamming)
  {
    femMembers = (clusters$raw_para == clusters.raw[fmod])
    consMembers = (clusters$hamming == clusters.hamming[cmod])
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    countbl[fmod, cmod] = sum(clusters$raw_para == clusters.raw[fmod] & clusters$hamming == clusters.hamming[cmod])
  }
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50

totals.raw = apply(countbl,1,sum)
totals.hamming = apply(countbl,2,sum)
pdf(file="hamming/consersus.pdf",h=8,w=14)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)

labeledHeatmap(Matrix = pTable,
yLabels = paste(" ", clusters.raw),
xLabels = paste(" ", clusters.hamming),
colorLabels = TRUE,
ySymbols = paste("ESCC_", clusters.raw, ": ", totals.raw, sep=""),
xSymbols = paste("EAC_", clusters.hamming, ": ", totals.hamming, sep=""),
textMatrix = countbl,
colors = greenWhiteRed(100)[50:100],
main = "Cell cluster overlap of default parameter and hamming",
cex.text = 1.0,
cex.lab = 1.0,
setStdMargins = FALSE)
dev.off()


pdf("other.vlnplot/KDR.pdf",h=4,w=6)
VlnPlot(data, features="KDR", pt.size=0)
dev.off()
