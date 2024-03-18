##R 
##single cell data of ESCC
##step2: EC sub-cluster
##Jiacheng Dai #daicy0424@gmail.com
##2023/04/26

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

#load
setwd("01ECatlas")
data = readRDS("01total_cluster/new0418/01ESCCdataset.211750.rds")
meta = read.csv("01total_cluster/new0418/metadata.211.csv", header=T, row.names=1)
meta.smc = meta[meta$celltype == "Smooth muscle cell",]
data.matrix = data@assays$RNA@counts[,rownames(meta.smc)]
# 28900 12252

data1 <- CreateSeuratObject(counts = data.matrix,
                        project = "SMC",
                        min.cells = 3,
                        min.features = 0)
options(future.globals.maxSize = 2000 * 1024^2)
data1 <- SCTransform(data1, verbose = T)

data1 = RunPCA(data1)
data1 = RunUMAP(data1, dims = 1:30)
data1 = FindNeighbors(data1, dims = 1:30) 
data1 = FindClusters(data1, resolution = 0.5) 

cluster.markers = FindAllMarkers(data1,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers, "08smc/SMCmarker.csv")

pdf("08smc/08umap.SCT05.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "SCT_snn_res.0.5", label=TRUE, pt.size = 0.5)
dev.off()
meta.smc$SMCcluster = data1@meta.data$seurat_clusters
write.csv(meta.smc, "08smc/metadata.smc.csv")

#
meta0 = meta.smc[meta.smc$SMCcluster == 0, ]; meta0$celltype = "SMC_1"
meta1 = meta.smc[meta.smc$SMCcluster == 1, ]; meta1$celltype = "SMC_2"
meta2 = meta.smc[meta.smc$SMCcluster == 2, ]; meta2$celltype = "Pericyte_1"
meta3 = meta.smc[meta.smc$SMCcluster == 3, ]; meta3$celltype = "SMC_3"
meta4 = meta.smc[meta.smc$SMCcluster == 4, ]; meta4$celltype = "SMC_4"
meta5 = meta.smc[meta.smc$SMCcluster == 5, ]; meta5$celltype = "SMC_5"
meta6 = meta.smc[meta.smc$SMCcluster == 6, ]; meta6$celltype = "SMC_6"
meta8 = meta.smc[meta.smc$SMCcluster == 8, ]; meta8$celltype = "SMC_8"
meta9 = meta.smc[meta.smc$SMCcluster == 9, ]; meta9$celltype = "Pericyte_2"
meta11 = meta.smc[meta.smc$SMCcluster == 11, ]; meta11$celltype = "SMC_9"
total = rbind(meta0, meta1, meta2, meta3, meta4, meta5, meta6, meta8, meta9, meta11)
data2 = data1[,rownames(total)]
data2@meta.data$celltype = total$celltype
data2@meta.data$tissue = total$tissue
#19108, 11435
saveRDS(data2, "08smc/08SMC.sct.rds")

pdf("08smc/08umap.cluster1.pdf",w=12,h=10)
DimPlot(data2, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()
pdf("08smc/08umap.celltype.pdf",w=6,h=5)
DimPlot(data2, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5)
dev.off()
pdf("08smc/08umap.tissuetype.pdf",w=12,h=10)
DimPlot(data2, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()

#Dotplot
cd_gene = c("CCL19","CCL21","CXCL12","CCL2","IGFBP5","FGF7","NDUFA4L2","COX4I2","B2M","CPE","GSN","IGFBP6","MT2A","STEAP4","FABP4","FABP5","IGFBP2","MFAP4","SFRP2","FBLN1","MMP2","HLA-DRA","HLA-DRB1","HSPA1A","HSPA1B","HSPA6","COL1A1","COL3A1","COL4A1","TIMP1",'RGS5',"SPARC","ISG15","HLA-A","HLA-B","HLA-C")
cd_gene = c("ACTA2","MYH11","MYL9","NR4A1","FOS","EGR1","ATF3","PLN","SORBS2","TAGLN","VIM","FLNA","CNN1","DES","ACTG2","RGS2","TPM1","GJA1","GREM1","SOD2","ID2","HOPX","FOXP1","CENPW")

data.usage = DotPlot(data2,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(1,4,5,6,7,11,2,0,3,8,9))
levels(data.usage$id) = c("Pericyte_1","Pericyte_2","Pericyte_3","Pericyte_4","Pericyte_5","Pericyte_6","T-Pericyte","SMC_1","SMC_2","SMC_3","T-SMC")

pdf("08dotplot.pericyte.pdf",h=5,w=14)
pdf("08dotplot.smc.pdf",h=5,w=10)
pdf("08dotplot.all.pdf",h=4.5,w=20)
data.usage %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=1)) +
ylab('') +
xlab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(0,3))
dev.off()

#
data=readRDS("08SMC.sct.rds")
meta = data3[[]]
meta0 = meta[meta$seurat_clusters == 0, ]; meta0$celltype = "SMC_1"
meta1 = meta[meta$seurat_clusters == 1, ]; meta1$celltype = "Pericyte_1"
meta2 = meta[meta$seurat_clusters == 2, ]; meta2$celltype = "T-Pericyte"
meta3 = meta[meta$seurat_clusters == 3, ]; meta3$celltype = "SMC_2"
meta4 = meta[meta$seurat_clusters == 4, ]; meta4$celltype = "Pericyte_2"
meta5 = meta[meta$seurat_clusters == 5, ]; meta5$celltype = "Pericyte_3"
meta6 = meta[meta$seurat_clusters == 6, ]; meta6$celltype = "Pericyte_4"
meta7 = meta[meta$seurat_clusters == 7, ]; meta7$celltype = "Pericyte_5"
meta8 = meta[meta$seurat_clusters == 8, ]; meta8$celltype = "SMC_3"
meta9 = meta[meta$seurat_clusters == 9, ]; meta9$celltype = "T-SMC"
meta11 = meta[meta$seurat_clusters == 11, ]; meta11$celltype = "Pericyte_6"
total = rbind(meta0, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta11)
data1 = data[,rownames(total)]
data1@meta.data$celltype = total$celltype

pdf("08umap.celltype.pdf",w=6.5,h=5)
DimPlot(data1, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
