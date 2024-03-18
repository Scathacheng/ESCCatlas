##R 
##single cell data of ESCC
##step5: compare lungEC data and esccEC data
##Jiacheng Dai #daicy0424@gmail.com
##2022/08/21

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(harmony)

setwd("01ECatlas/05lungEC")

data = read.csv("Data.csv", header=T, row.names=1)
data1 = data.frame(data)
lungECdata = as(as.matrix(data1),"dgCMatrix") #14917 21252

escc = readRDS("../02EC_cluster/new0423/02ECdataset.16580.rds")
data.matrix = escc@assays$SCT@data
genes = intersect(rownames(lungECdata),rownames(data.matrix)) #14315

metadata_lung = read.csv("Metadata.csv", header=T, row.names=1)
metadata_escc = read.csv("total.16580.csv", header=T, row.names=1)
meta1 = metadata_lung[,c("Tumor.type","NEC.TEC","Study","Cluster")]
meta2 = metadata_escc[,c("orig.ident","tissue","orig.ident","celltype")]
colnames(meta1) = colnames(meta2) = c("Tumortype","Origin","Study","Cluster")
meta2$Tumortype = "squamous"
meta2$Study = "esophagusEC"
combineMeta = rbind(meta1, meta2)
write.csv(combineMeta,"05combineMeta.csv")

lungECdata1 = lungECdata[genes,]
esccECdata1 = data.matrix[genes,rownames(metadata_ec)]
lungECdata2 = scale(lungECdata1)
esccECdata2 = scale(esccECdata1)
data = cbind(lungECdata2, esccECdata2)

###################################################################################################################
#lungEC solo
seuratObject = CreateSeuratObject(counts=lungECdata)
seuratObject <- FindVariableFeatures(object = seuratObject, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE, x.low.cutoff = 0.00125, x.high.cutoff = 8, y.cutoff = 0.25, do.cpp = TRUE)
metadata_lung = read.csv("Metadata.csv", header=T, row.names=1)
seuratObject@meta.data[4:8] = metadata_lung[1:5]

all.genes = rownames(seuratObject)
seuratObject = ScaleData(seuratObject, features = all.genes)
seuratObject = RunPCA(seuratObject, features=VariableFeatures(object=seuratObject),npcs=15)

seuratObject = RunUMAP(seuratObject, dims=1:15)
seuratObject = FindNeighbors(seuratObject, dims=1:15)
seuratObject = FindClusters(seuratObject, reduction.type = "custom", resolution=1, k.param = 10, dims = 1:15, random.seed =1234)

pdf("UMAP_lungEC_celltype.pdf", h=5,w=8)
DimPlot(seuratObject, reduction="umap", label=TRUE, label.size=4, group.by="Cluster", pt.size=.1)
dev.off()
pdf("UMAP_lungEC_Cluster.1.pdf", h=5,w=7)
DimPlot(seuratObject, reduction="umap", label=TRUE, pt.size=.1)
dev.off()

saveRDS(seuratObject, "05lungEC_res.1.rds")

pdf("01lungEC_solo/FeaturePlot.PROX1.pdf",h=5,w=7)
FeaturePlot(seuratObject,features=c("PROX1"))
dev.off()
#####################################################################################################################


combineObject = CreateSeuratObject(counts = data)
#saveRDS(combineObject, "05combinedata.zscore.rds")
combineObject@meta.data[,4:7] = combineMeta[,1:4]
#combineObject = NormalizeData(combineObject, normalization.method="LogNormalize", scale.factor=10000)
combineObject = FindVariableFeatures(combineObject, selection.method = "vst", nfeatures=2000)

all.genes = rownames(combineObject)
combineObject = ScaleData(combineObject, features = all.genes)
combineObject = RunPCA(combineObject, features = VariableFeatures(object = combineObject), npcs = 15)

#plot before harmony
p1 = DimPlot(object = combineObject, reduction = "pca", pt.size = .1, group.by = "Study")
p2 = VlnPlot(object = combineObject, features = "PC_1", pt.size = 0, group.by = "Study")
pdf("02harmony/plot_bfHarmony.pdf",h=5,w=10)
p1+p2
dev.off()

#Run harmony
combineObject = RunHarmony(combineObject, "Study", plot_comvergence=TRUE)
harmony_embeddings = Embeddings(combineObject, 'harmony')
harmony_embeddings[1:5,1:5]

p1 = DimPlot(object = combineObject, reduction = "harmony", pt.size = .1, group.by = "Study")
p2 = VlnPlot(object = combineObject, features = "harmony_1", pt.size = 0, group.by = "Study")
pdf("02harmony/plot_aftHarmony.pdf",h=5,w=10)
p1+p2
dev.off()

#Downstream analysis
combineObject = RunUMAP(combineObject, reduction="harmony", dims=1:15)
combineObject = FindNeighbors(combineObject, dims=1:15)
combineObject = FindClusters(combineObject, reduction.type = "custom", resolution=2, k.param = 10, dims = 1:15, random.seed =1234)
combineObject.markers = FindAllMarkers(combineObject, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

pdf("UMAP_Study.pdf", h=5,w=7)
DimPlot(combineObject, reduction="umap", group.by="Study", pt.size=.1)
dev.off()
pdf("UMAP_Origin.pdf", h=5,w=7)
DimPlot(combineObject, reduction="umap", group.by="Origin", pt.size=.1)
dev.off()
pdf("UMAP_Type.pdf", h=5,w=7)
DimPlot(combineObject, reduction="umap", group.by="Tumortype", pt.size=.1)
dev.off()
pdf("UMAP_SeuratCluster_res2.pdf", h=5,w=7)
DimPlot(combineObject, reduction="umap", label=TRUE, pt.size=.1, repel=TRUE)
dev.off()

saveRDS(combineObject, "05combinedata.res2.rds")
write.csv(combineObject@meta.data, "05combinedata.metadata.csv")

pdf("UMAP_lungOrigin.pdf",h=5,w=6)
DimPlot(data2,reduction="umap", label=FALSE, pt.size=1, group.by="Origin", cols=c('TEC'='#F8766D','NEC'='#00BA38'))
dev.off()
pdf("UMAP_esoOrigin.pdf",h=5,w=6)
DimPlot(data1,reduction="umap", label=FALSE, pt.size=1, group.by="Origin", cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()
pdf("UMAP_cellident2.pdf", h=5,w=8)
DimPlot(data, reduction="umap", label=TRUE, label.size=4, repel=TRUE, group.by="cell.ident", pt.size=.1)
dev.off()


#
combineObject = readRDS("05combinedata.res2.rds")
meta = combineObject[[]]
meta1 = meta[meta$seurat_clusters %in% c(1,8),]; meta1$cell.ident = "artery"
meta2 = meta[meta$seurat_clusters %in% c(2,26,22),]; meta2$cell.ident = "capillary"
meta3 = meta[meta$seurat_clusters %in% c(3,12,20,24,9),]; meta3$cell.ident = "capillary Lung specific"
meta4 = meta[meta$seurat_clusters %in% c(21,23,33,11,6,28),]; meta4$cell.ident = "lymphatic EC"
meta5 = meta[meta$seurat_clusters %in% c(7,16,13,17,31,25,14,19,27),]; meta5$cell.ident = "vein"
meta6 = meta[meta$seurat_clusters %in% c(4,15),]; meta6$cell.ident = "stalk cell"
meta7 = meta[meta$seurat_clusters %in% c(10,18,29,32),]; meta7$cell.ident = "tip cell"
meta8 = meta[meta$seurat_clusters %in% c(5,0,30),]; meta8$cell.ident = "immunomodulatory EC"
total = rbind(meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8)
total = total[rownames(combineObject@meta.data),]
write.csv(total, "metadata.40389.csv")

combineObject@meta.data$cell.ident = total$cell.ident
combineObject@meta.data$Cluster = total$Cluster

pdf("UMAP_cellident.pdf", h=10,w=15)
DimPlot(combineObject, reduction="umap", label=TRUE, label.size=4, repel=TRUE, group.by="cell.ident", pt.size=.1)
dev.off()
pdf("UMAP_celltype.pdf", h=10,w=18)
DimPlot(combineObject, reduction="umap", label=TRUE, label.size=4, repel=TRUE, group.by="Cluster", pt.size=.1)
dev.off()

#Feature plot
pdf("03FeaturePlot/COL4A1.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("COL4A1"))
dev.off()
pdf("03FeaturePlot/INSR.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("INSR"))
dev.off()
pdf("03FeaturePlot/GJA5.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("GJA5"))
dev.off()
pdf("03FeaturePlot/ACKR1.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("ACKR1"))
dev.off()
pdf("03FeaturePlot/CA4.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("CA4"))
dev.off()
pdf("03FeaturePlot/PROX1.pdf",w=6,h=5)
FeaturePlot(combineObject,features=c("PROX1"))
dev.off()

#AverageExpression
meta = read.csv("marker_res2.csv",header=T)
average = AverageExpression(combineObject, group.by = "Cluster", feature=unique(meta$gene))
test = average$RNA
clusterTree = hclust(dist(t(test),method="canberra"),method="complete")
pdf("clusterTree.pdf")
plot(clusterTree)
dev.off()

average = AverageExpression(combineObject, group.by = "Cluster")
test = average$RNA
test1 = log(test+1,base=2)
cor.data = corr.test(test1[,c(1,2,6,7,8,9,10,16,18,19,20,22,24)], test1[,c(3,4,5,11,12,13,14,15,17,21,23,25,26,27,28,29,30)], method="pearson")
data=t(cor.data$r)
p=t(cor.data$p.adj)

pdf("03FeaturePlot/CorHeatmap.pdf",h=4.5,w=7)
pheatmap(data, display_numbers = matrix(ifelse(data > 0.85, "*", ''), nrow(data)))
dev.off()

#
options(stringsAsFactors=TRUE)
meta = read.csv("05combinedata.metadata.csv",header=T,row.names=1)

proplist=matrix(NA,nrow=18,ncol=2)
j=1
for (i in 0:17) {
    meta1 = meta[meta$RNA_snn_res.1==i,]
    proplist[j,] = table(meta1$Study)
    j = j+1
}

table(meta1 = meta[meta$RNA_snn_res.1==3,]$Cluster)
table(meta[meta$Cluster=="CCL2+ vein",]$RNA_snn_res.1)

#dotplot
library(Seurat)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(patchwork)
combineObject = readRDS("05combinedata.rds")
#all, h=7,w=8
cd_gene<-c("PECAM1","ACKR1","VWF","SOX17","EFNB2","GJA5","GJA4","CAV1","CA4","LYVE1","PROX1","ENG","HSPG2","INSR","PLVAP","APLNR","CXCR4","PGF","LXN","COL4A1","COL4A2")
#tip_cell_combine, h=6,w=18
cd_gene<-c("ADM","ANGPT2","APLN","CXCR4","ESM1","PXDN","LXN","ANGPTL2","INSR","FILIP1","ACTN4","MYH9","MYL9","MYO1B","FSCN1","MARCKS","MARCKSL1","CD93","LGALS1","ITGB1","MCAM","SPARC","SPARCL1","LAMA4","LAMC1","COL4A1","COL4A2","SERPINH1","LOX","CTHRC1","ENG","PLVAP","HSPG2","APLNR","VWA1","RBP7","MMP2","KDR","FLT1","TIE1")

cd_gene<-c("GJA5","EFNB2","SOX17","CXCL12","DKK2","FBLN2","MGP","SRGN","IL33","FN1")
cd_gene<-c("ACKR1","SELP","VCAM1","VWF","CLU","CPE","CCL15","CCL23","NAMPT","NNMT")
cd_gene<-c("FCN3","SLC6A4","TMEM100","VIPR1","CAV1","NOSTRIN","CA4","FABP5")
cd_gene<-c("LYVE1","CCL21","TFF3","MMRN1","ADIRF","AKAP12","CD9","ANXA2","SNCG","GYPC","S100A10","PPFIBP1","EFEMP1")
#activePCV, h=6,w=7
cd_gene<-c("ACKR1","SELP","VCAM1","POSTN","IGFBP7","CCL14","PRCP","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DRA","HLA-DRB1")

data.usage = DotPlot(combineObject,features=cd_gene)$data
data.usage$id2 = factor(data.usage$id, levels=c(5,16,0,10,12,13,14,17,2,8,11,15,4,9,6,7,1,3))
levels(data.usage$id2) = c("artery1 (C5_ESCC)","artery2 (C16_NSCLC)","vein1 (C0_ESCC)","vein2 (C10_ESCC)","vein3 (C12_ESCC)","vein4 (C13_NSCLC)","vein5 (C14_ESCC)","vein6 (C17_ESCC)","capillary1 (C2_NSCLC)","capillary2 (C8_NSCLC)","capillary3 (C11_NSCLC)","capillary4 (C15_ESCC)","TEC (C4_ESCC)","tip cell (C9_NSCLC)","LEC1 (C6_ESCC)","LEC2 (C7_NSCLC)","Qsc EC (C1_ESCC)","RGCC+ EC (C3_ESCC)")

pdf("03FeaturePlot1019/dotplot_cap.pdf",h=7,w=8)
ggplot(data.usage, aes(x=features.plot, y=id2, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
ylab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(0,3))
dev.off()

#Sankey plot
library(tidyverse)
library(viridis)
library(patchwork)
library(networkD3)
data=read.csv("04SankeyPlot/data.csv",header=T)
nodes=read.csv("04SankeyPlot/node.csv",header=T)
sankeyNetwork(Links = data, Nodes=nodes, Source = "IDsource", Target="IDtarget", Value="value1", NodeID="names",sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

#heatmap
library(Seurat)
library(psych)
library(pheatmap)
combineObject = readRDS("05combinedata.rds")
combineObject@active.ident = factor(combineObject$Cluster)
levels(combineObject@active.ident)[24] = "KRT+ EC"
meta = read.csv("05markers.1.csv",header=T)
average = AverageExpression(combineObject,feature=unique(meta$gene))
cor.data = corr.test(average$RNA[,c(1,2,6,7,8,9,10,16,18,19,20,22,24)], average$RNA[,c(3,4,5,11,12,13,14,15,17,21,23,25,26,27,28,29,30)], method="pearson")
data=t(cor.data$r)
p=t(cor.data$p.adj)

pdf("03FeaturePlot1019/CorHeatmap1.pdf",h=4.5,w=7)
pheatmap(data, display_numbers = matrix(ifelse(data > 0.85, "*", ''), nrow(data)))
dev.off()

#SCENIC
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(AUCell)

dir.create("int")
meta = combineObject[[]]
meta = meta[meta$Cluster %in% c("activated pcv","arteries","capillary scavenging","capillary intermediate","tip cell","immature"),]

cellinfo = meta[,c("Cluster","cell.ident")]
saveRDS(cellinfo, file="int/cellInfo.rds")

exprMat = as.matrix(combineObject@assays$RNA@counts[,rownames(cellinfo)])
mydbDIR = "./cisTarget"
mydbs = c("hg19-tss-centered-10kb-7species.mc9nr.feather","hg19-500bp-upstream-7species.mc9nr.feather")
names(mydbs) = c("10kb","500bp")
scenicOptions = initializeScenic(org="hgnc",
                                nCores=16,
                                dbDir=mydbDIR,
                                dbs=mydbs,
                                datasetTitle="lungEC")
saveRDS(scenicOptions, "int/scenicOp-lungEC.rds")

genesKept = geneFiltering(exprMat, scenicOptions, minCountsPerGene = 3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+5) 
runGenie3(exprMat_filtered_log, scenicOptions)

saveRDS(cellinfo, file="int/cellInfo.rds")
saveRDS(scenicOptions, "int/scenicOp-lungEC.rds")

exprMat_log <- log2(exprMat+6)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

saveRDS(cellinfo, file="int/cellInfo.rds")
saveRDS(scenicOptions, "int/scenicOp-lungEC.rds")

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(aucell_regulonAUC),cellAnnotation=cellinfo[colnames(aucell_regulonAUC),"Cluster"])
rssPlot <- plotRSS(rss)
rssPlot$plot
pdf("rssplot.pdf")
rssPlot$plot
dev.off()

regulonAUC=loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC=regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType = sapply(split(rownames(cellinfo),cellinfo$Cluster),function(cells)rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("ComplexHeatmap.pdf")
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()
