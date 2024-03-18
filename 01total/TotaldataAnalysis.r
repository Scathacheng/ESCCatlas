#

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(harmony)

#load data
setwd("new0418")
meta = read.csv("metadata.211.csv",header=T,row.names=1)
data = readRDS("01ESCCdataset.211750.rds")


#filtering sample
#meta1 = meta[meta$Patient_ID %in% c("p026","p036","p037","p041","p045","p047","p048","p049","p050","p054","p055"),] 
meta1 = meta[meta$celltype != "Fibroblast.independent",]
data1 = data[,rownames(meta1)]
#138356


data1@meta.data$batch = meta1$Sample_batch
data1@meta.data$patient = meta1$Patient_ID

#harmony
data = readRDS("01Alldata.rds")
data = RunHarmony(data, "batch", assay.use = "SCT", project.dim=FALSE)

p1 = DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "batch")
p2 = VlnPlot(object = data, features = "PC_1", pt.size = 0, group.by = "batch")

p3 = DimPlot(object = data, reduction = "harmony", pt.size = .1, group.by = "batch")
p4 = VlnPlot(object = data, features = "harmony_1", pt.size = 0, group.by = "batch")

pdf("batch-pc.pdf",h=5,w=10)
p1+p2
dev.off()

pdf("batch-harmony.pdf",h=5,w=10)
p3+p4
dev.off()

data = RunUMAP(data, reduction="harmony", dims=1:15)
data = FindNeighbors(data, dims=1:15)
data = FindClusters(data, reduction.type = "custom", resolution=1, k.param = 10, dims = 1:15, random.seed = 888)
saveRDS(data, "01Alldata.harmony.rds")
markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "marker100.csv")
saveRDS(data, "01Alldata.harmony.rds")

DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)

meta = data[[]]
meta$major = NA
meta[meta$seurat_clusters %in% c(2,17),]$major = "B cell"
meta[meta$seurat_clusters %in% c(0,1,7,9,23),]$major = "Fibroblast"
meta[meta$seurat_clusters %in% c(3,20,27),]$major = "Endothelial cell"
meta[meta$seurat_clusters==29,]$major = "Lymphatic"
meta[meta$seurat_clusters %in% c(4,5,6,12,13,22),]$major = "T cell"
meta[meta$seurat_clusters %in% c(8,11,19,21,37,38),]$major = "Epithelial cell"
meta[meta$seurat_clusters %in% c(14,15,25,26),]$major = "Myeloid cell"
meta[meta$seurat_clusters %in% c(10,16),]$major = "Smooth muscle cell"
meta[meta$seurat_clusters== 18,]$major = "Mast cell"
meta[meta$seurat_clusters %in% c(24,28,30,31,32,33,34,35,36,39),]$major = "Plasma cell"
write.csv(meta, "meta_allsample.csv")
data@meta.data$celltype = meta$major

pdf("All_tissue.pdf",w=8,h=6)
DimPlot(data, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.001, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()
pdf("All_celltype.pdf",w=9,h=6)
DimPlot(data, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("All_patient.pdf",w=7,h=6)
DimPlot(data, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.1)
dev.off()
pdf("All_batch.pdf",w=7,h=6)
DimPlot(data, reduction = "umap", group.by = "batch", label=FALSE, pt.size = 0.1)
dev.off()
pdf("All_seuratcluster.pdf",w=8,h=6)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()

##downsampling
setwd("ds25")
data = readRDS("../01Alldata.harmony.rds")
# 25%
set.seed(1234)
dws.25 = data[,sample(colnames(data),size = 51595, replace=F)]
dws.25 = RunUMAP(dws.25, reduction = "harmony", dims=1:15)
dws.25 = FindNeighbors(dws.25, dims=1:15)
dws.25 = FindClusters(dws.25, reduction.type = "custom", resolution=1, k.param = 10, dims = 1:15, random.seed =888)
markers = FindAllMarkers(dws.25, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(dws.25, "sample25.rds")
write.csv(markers, "marker25.csv")

meta = dws.25[[]]
meta$major = NA
meta[meta$seurat_clusters== 0,]$major = "B cell"
meta[meta$seurat_clusters %in% c(1,4,6,8),]$major = "Fibroblast"
meta[meta$seurat_clusters %in% c(5,12),]$major = "Endothelial cell"
meta[meta$seurat_clusters== 27,]$major = "Lymphatic"
meta[meta$seurat_clusters %in% c(3,7,9,11,15,18),]$major = "T cell"
meta[meta$seurat_clusters %in% c(2,17,21,23,30),]$major = "Epithelial cell"
meta[meta$seurat_clusters %in% c(10,14,26),]$major = "Myeloid cell"
meta[meta$seurat_clusters %in% c(13,20),]$major = "Smooth muscle cell"
meta[meta$seurat_clusters== 16,]$major = "Mast cell"
meta[meta$seurat_clusters %in% c(19,22,24,25,28,29),]$major = "Plasma cell"

write.csv(meta, "meta_ds25.csv")
dws.25@meta.data$celltype = meta$major

pdf("ds25_tissue.pdf",w=8,h=6)
DimPlot(dws.25, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.001, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()
pdf("ds25_celltype.pdf",w=9,h=6)
DimPlot(dws.25, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("ds25_patient.pdf",w=7,h=6)
DimPlot(dws.25, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds25_batch.pdf",w=7,h=6)
DimPlot(dws.25, reduction = "umap", group.by = "batch", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds25_seuratcluster.pdf",w=8,h=6)
DimPlot(dws.25, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()

# 50%
setwd("ds50")
data = readRDS("../01Alldata.harmony.rds")
set.seed(2345)
dws.50 = data[,sample(colnames(data),size = 103190, replace=F)]
dws.50 = RunUMAP(dws.50, reduction = "harmony", dims=1:15)
dws.50 = FindNeighbors(dws.50, dims=1:15)
dws.50 = FindClusters(dws.50, reduction.type = "custom", resolution=1, k.param = 10, dims = 1:15, random.seed = 888)
markers = FindAllMarkers(dws.50, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(dws.50, "sample50.rds")
write.csv(markers, "marker50.csv")

meta = dws.50[[]]
meta$major = NA
meta[meta$seurat_clusters== 0,]$major = "B cell"
meta[meta$seurat_clusters %in% c(1,2,8,9,25),]$major = "Fibroblast"
meta[meta$seurat_clusters %in% c(3,10),]$major = "Endothelial cell"
meta[meta$seurat_clusters== 28,]$major = "Lymphatic"
meta[meta$seurat_clusters %in% c(4,6,7,11,13,21),]$major = "T cell"
meta[meta$seurat_clusters %in% c(5,12,20,24),]$major = "Epithelial cell"
meta[meta$seurat_clusters %in% c(15,16,17),]$major = "Myeloid cell"
meta[meta$seurat_clusters %in% c(14,18),]$major = "Smooth muscle cell"
meta[meta$seurat_clusters== 19,]$major = "Mast cell"
meta[meta$seurat_clusters %in% c(22,23,26,27,29,30,31,32),]$major = "Plasma cell"
write.csv(meta, "meta_ds50.csv")
dws.50@meta.data$celltype = meta$major

pdf("ds50_tissue.pdf",w=8,h=6)
DimPlot(dws.50, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.001, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()
pdf("ds50_celltype.pdf",w=9,h=6)
DimPlot(dws.50, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("ds50_patient.pdf",w=7,h=6)
DimPlot(dws.50, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds50_batch.pdf",w=7,h=6)
DimPlot(dws.50, reduction = "umap", group.by = "batch", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds50_seuratcluster.pdf",w=8,h=6)
DimPlot(dws.50, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()

# 75%
setwd("ds75")
data = readRDS("../01Alldata.harmony.rds")
set.seed(3456)
dws.75 = data[,sample(colnames(data),size = 154785, replace=F)]
dws.75 = RunUMAP(dws.75, reduction = "harmony", dims=1:15)
dws.75 = FindNeighbors(dws.75, dims=1:15)
dws.75 = FindClusters(dws.75, reduction.type = "custom", resolution=1, k.param = 10, dims = 1:15, random.seed = 888)
markers = FindAllMarkers(dws.75, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(dws.75, "sample75.rds")
write.csv(markers, "marker75.csv")

meta = dws.75[[]]
meta$major = NA
meta[meta$seurat_clusters== 0,]$major = "B cell"
meta[meta$seurat_clusters %in% c(1,4,7,12,20,21),]$major = "Fibroblast"
meta[meta$seurat_clusters %in% c(2,18,19),]$major = "Endothelial cell"
meta[meta$seurat_clusters== 32,]$major = "Lymphatic"
meta[meta$seurat_clusters %in% c(3,5,6,8,10,29),]$major = "T cell"
meta[meta$seurat_clusters %in% c(9,11,22,23,27),]$major = "Epithelial cell"
meta[meta$seurat_clusters %in% c(15,16,24,30),]$major = "Myeloid cell"
meta[meta$seurat_clusters %in% c(13,14),]$major = "Smooth muscle cell"
meta[meta$seurat_clusters== 17,]$major = "Mast cell"
meta[meta$seurat_clusters %in% c(25,26,28,31,33,34,35,36,37,38),]$major = "Plasma cell"
write.csv(meta, "meta_ds75.csv")
dws.75@meta.data$celltype = meta$major

pdf("ds75_tissue.pdf",w=8,h=6)
DimPlot(dws.75, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.1, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()
pdf("ds75_celltype.pdf",w=9,h=6)
DimPlot(dws.75, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("ds75_patient.pdf",w=7,h=6)
DimPlot(dws.75, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds75_batch.pdf",w=7,h=6)
DimPlot(dws.75, reduction = "umap", group.by = "batch", label=FALSE, pt.size = 0.1)
dev.off()
pdf("ds75_seuratcluster.pdf",w=8,h=6)
DimPlot(dws.75, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()


##overlap
data1 = read.csv("meta_allsample.csv",header=T,row.names=1)
data1$cellID = rownames(data1)
ct.ref = unique(data1$major)
me.ref = data1[,c("cellID","major")]

data2 = read.csv("meta_ds75.csv",header=T,row.names=1) 
data2$cellID = rownames(data2)
ct.test = unique(data2$major)
me.test = data2[,c("cellID","major")]

pTable = matrix(0, nrow=10, ncol=10)
countbl = matrix(0, nrow=10, ncol=10)

for (fmod in 1:10)
    for (cmod in 1:10)  
    {
        femMembers = me.ref[me.ref$major == ct.ref[fmod],]$cellID
        consMembers = me.test[me.test$major == ct.test[cmod],]$cellID
        a1 = length(intersect(femMembers,consMembers))
        b1 = length(femMembers) - a1
        c1 = length(consMembers) - a1
        d1 = 206380-a1-b1-c1
        pTable[fmod, cmod] = -log10(fisher.test(rbind(c(a1,b1),c(c1,d1)), alternative = "greater")$p.value)
        countbl[fmod, cmod] = a1
    }
#num=1.3*max(pTable[is.finite(pTable)])
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]) 
rownames(pTable)  = rownames(countbl) = ct.ref
colnames(pTable) = colnames(countbl) = ct.test

totals.ref = apply(countbl,1,sum)
totals.test = apply(countbl,2,sum)

par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(10, 13, 2.7, 1)+0.3)
labeledHeatmap(Matrix = pTable,
yLabels = paste("100%_", ct.ref, ": ", totals.ref, sep=""),
xLabels = paste("75%_", ct.test, ": ", totals.test, sep=""),
colorLabels = TRUE,
#ySymbols = paste("100%_", ct.ref, ": ", totals.ref, sep=""),
#xSymbols = paste("25%_", ct.test, ": ", totals.test, sep=""),
textMatrix = countbl,
colors = blueWhiteRed(100)[50:100],
main = "Down Sampling Overlap",
cex.text = 1.0,
cex.lab = 1.0,
setStdMargins = FALSE)


#feature plot
pdf("featureplot/endo.pdf",h=4,w=4.5)
FeaturePlot(data,features="CDH5")
dev.off()

pdf("featureplot/smc.pdf",h=4,w=4.5)
FeaturePlot(data,features="ACTA2")
dev.off()

pdf("featureplot/mast.pdf",h=4,w=4.5)
FeaturePlot(data,features="TPSB2")
dev.off()

pdf("featureplot/fibro.pdf",h=4,w=4.5)
FeaturePlot(data,features="DCN")
dev.off()

pdf("featureplot/epithelial.pdf",h=4,w=4.5)
FeaturePlot(data,features="KRT6A")
dev.off()

pdf("featureplot/tcell.pdf",h=4,w=4.5)
FeaturePlot(data,features="CD3D")
dev.off()

pdf("featureplot/bcell.pdf",h=4,w=4.5)
FeaturePlot(data,features="CD79A")
dev.off()

pdf("featureplot/myeloid1.pdf",h=4,w=4.5)
FeaturePlot(data,features="C1QA")
dev.off()

pdf("featureplot/myeloid2.pdf",h=4,w=4.5)
FeaturePlot(data,features="CXCL8")
dev.off()

pdf("featureplot/lec.pdf",h=4,w=4.5)
FeaturePlot(data,features="CCL21")
dev.off()

#re-clustering
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
setwd("03GeosAltas")
meta = read.csv("new1228/ds100/meta_allsample.csv",header=T,row.names=1)
data = readRDS("new1228/ds100/01Alldata.harmony.rds")
data@meta.data$celltype = meta$major

meta2 = meta[meta$major == "Epithelial cell",]
data1 = data[,rownames(meta2)]
saveRDS(data1, "02epithelial/02data.23652.rds")
write.csv(meta2, "02epithelial/metadata.23652.csv")

meta2 = meta[meta$major == "T cell",]
data1 = data[,rownames(meta2)]
saveRDS(data1, "03tcell/03data.47724.rds")
write.csv(meta2, "03tcell/metadata.47724.csv")

meta2 = meta[meta$major %in% c("B cell","Plasma cell"),]
data1 = data[,rownames(meta2)]
saveRDS(data1, "04bcell/04data.33202.rds")
write.csv(meta2, "04bcell/metadata.33202.csv")

meta2 = meta[meta$major == "Fibroblast",]
data1 = data[,rownames(meta2)]
saveRDS(data1, "05fibroblast/05data.47920.rds")
write.csv(meta2, "05fibroblast/metadata.47920.csv")

meta2 = meta[meta$major == "Myeloid cell",]
data1 = data[,rownames(meta2)]
saveRDS(data1, "06myeloidcell/06data.16636.rds")
write.csv(meta2, "06myeloidcell/metadata.16636.csv")

meta2 = meta[meta$major == "Endothelial cell",]
data1 = data[,rownames(meta2)]
saveRDS(data1, "07endothelial/07data.19288.rds")
write.csv(meta2, "07endothelial/metadata.19288.csv")
