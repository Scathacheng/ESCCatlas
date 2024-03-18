##R 
##single cell data of ESCC
##step2: fibroblast
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
meta.fb = meta[meta$celltype == "Fibroblast",]
data.matrix = data@assays$RNA@counts[,rownames(meta.fb)]
#28900 46870

data1 <- CreateSeuratObject(counts = data.matrix,
                        project = "fibroblast",
                        min.cells = 3,
                        min.features = 0)
options(future.globals.maxSize = 2000 * 1024^2)
data1 <- SCTransform(data1, verbose = T)

#clustering
data1 = RunPCA(data1)
data1 = RunUMAP(data1, dims = 1:30)
data1 = FindNeighbors(data1, dims = 1:30) 
data1 = FindClusters(data1, resolution = 0.3) 
saveRDS(data1, "08fibroblast/08FB.sct.rds")

#Find cluster biomarkers
cluster.markers = FindAllMarkers(data1,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers, "08fibroblast/FBmarker.res03.csv")

#plot 
pdf("08fibroblast/08umap.SCT1.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "SCT_snn_res.1", label=TRUE, pt.size = 0.5)
dev.off()
meta.fb$FBcluster = data1@meta.data$seurat_clusters
write.csv(meta.fb, "08fibroblast/metadata.fb.csv")

#
meta0 = meta.fb[meta.fb$seurat_clusters == 0, ]; meta0$celltype = "FB_1"
meta1 = meta.fb[meta.fb$seurat_clusters == 1, ]; meta1$celltype = "FB_2"
meta2 = meta.fb[meta.fb$seurat_clusters == 2, ]; meta2$celltype = "FB_3"
meta3 = meta.fb[meta.fb$seurat_clusters == 3, ]; meta3$celltype = "FB_4"
meta4 = meta.fb[meta.fb$seurat_clusters == 4, ]; meta4$celltype = "CAF_1"
meta5 = meta.fb[meta.fb$seurat_clusters == 5, ]; meta5$celltype = "CAF_2"
meta6 = meta.fb[meta.fb$seurat_clusters == 6, ]; meta6$celltype = "FB_5"
meta7 = meta.fb[meta.fb$seurat_clusters == 7, ]; meta7$celltype = "CAF_3"
meta8 = meta.fb[meta.fb$seurat_clusters == 8, ]; meta8$celltype = "FB_6"
meta9 = meta.fb[meta.fb$seurat_clusters == 9, ]; meta9$celltype = "FB_7"

total = rbind(meta0, meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9)
data2 = data1[,rownames(total)]
data2@meta.data$celltype = total$celltype
data2@meta.data$tissue = total$tissue
saveRDS(data2, "08fibroblast/08FB.sct.rds")

#plot 
pdf("08fibroblast/08umap.cluster03.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()
pdf("08fibroblast/08umap.celltype.pdf",w=12,h=10)
DimPlot(data1, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5)
dev.off()
pdf("08fibroblast/08umap.tissuetype.pdf",w=6,h=5)
DimPlot(data1, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#619CFF','Distant'='#00BA38'))
dev.off()





#DotPlot
library(Seurat)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggplot2)
library(patchwork)
cd_gene = c("ACTA2","S100A4","VIM","FAP","PDGFRA","PDGFRB")
data.usage = DotPlot(data,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(4,5,7,0,1,2,3,6,8,9))
levels(data.usage$id) = c("CAF_1","CAF_2","CAF_3","Fib_1","Fib_2","Fib_3","Fib_4","Fib_5","Fib_6","Fib_7")
pdf("08fibroblast/dotplot_activation.pdf",h=6,w=8)
data.usage %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=1)) +
ylab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(0,3))
dev.off()

#CAF
cd_gene = c("ACTA2","FAP","PDGFRB","POSTN","COL1A1","COL3A1","TPM2","TAGLN","MMP11","MMP1","MMP3","MMP13","CXCL8","SPARC","FN1","S100A9","SRGN","FABP5","C1QA","C1QB","HLA-DRA","HLA-DPA1","HLA-DPB1","CXCR4","TOP2A","MKI67","CENPF","UBE2C","TIMP1","IGFBP7","LXN","TGM2")
data.usage = DotPlot(data,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(4,5,7,0,1,2,3,6,8,9))
levels(data.usage$id) = c("CAF_1","CAF_2","CAF_3","Fib_1","Fib_2","Fib_3","Fib_4","Fib_5","Fib_6","Fib_7")
pdf("08fibroblast/dotplot_CAF.pdf",h=4.5,w=11)
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

#fibroblast
cd_gene = c("PTGDS","CCL13","CCL11","IGFBP4","IGFBP5","SFRP1","S100A10","MFAP5","SLPI","CCL2","CXCL2","JUNB","SOD2","PLA2G2A","RARRES1","TIMP1","SCN7A","SPARCL1","TCF21","PDGFRA","HSPA1A","HSPA1B","HSP90AA1","CXCL14","TGFBI","PDGFRL","CAV1","VIT","CYP1B1","SFRP4")
data.usage = DotPlot(data,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(4,5,7,0,1,2,3,6,8,9))
levels(data.usage$id) = c("CAF_1","CAF_2","CAF_3","Fib_1","Fib_2","Fib_3","Fib_4","Fib_5","Fib_6","Fib_7")
pdf("08fibroblast/dotplot_Fib.pdf",h=4.5,w=12)
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
