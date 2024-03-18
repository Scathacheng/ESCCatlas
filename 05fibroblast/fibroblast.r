##R 
##single cell data of ESCC
##step2: fibroblast
##Jiacheng Dai #daicy0424@gmail.com
##2023/04/26

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ROGUE)

#load data
data = readRDS("05data.47920.rds")
meta = read.csv("metadata.47920.csv",header=T,row.names=1)
#select resolution using ROGUE
expr = data@assays$RNA@counts
expr = as.matrix(expr)

data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
list = seq(from = 0.1, to = 3, by = 0.1) 
all.rogue = c()
for (i in 1:30){
    data = FindClusters(object = data, resolution = list[i])
    meta1 = data[[]]
    meta$ct = meta1[,9]
    rogue.res = rogue(expr, 
            labels = meta$ct, 
            samples = meta$patient, 
            platform = "UMI", 
            remove.outlier.n = 5,
            filter = T,
            min.cell.n = 300)
    av.rogue <- c()
    for (j in 1:ncol(rogue.res)) {
    tmp.r <- rogue.res[,j]
    tmp.r <- tmp.r[!is.na(tmp.r)]
    av.rogue[j] <- mean(tmp.r)
    }
    all.rogue[i] = mean(na.omit(av.rogue))
}
all.rogue

#res=1.3
data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(object=data, resolution=1.3, reduction.type = "custom", k.param = 10, dims = 1:15, random.seed =888)
saveRDS(data, "05data.res1.3.rds")


#Find cluster biomarkers
cluster.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers, "05markers.csv")

meta$fb_cluster = data@meta.data$seurat_clusters
write.csv(meta, "metadata.res1.3.csv")

#
meta$fb_celltype="NA"
meta[meta$fb_cluster==0,]$fb_celltype = "Fib.C0.MGP"
meta[meta$fb_cluster==1,]$fb_celltype = "Fib.C1.SLPI"
meta[meta$fb_cluster==2,]$fb_celltype = "Fib.C2.PTX3"
meta[meta$fb_cluster==3,]$fb_celltype = "Fib.C3.PTGDS"
meta[meta$fb_cluster==4,]$fb_celltype = "Fib.C4.PTGDS"
meta[meta$fb_cluster==5,]$fb_celltype = "Fib.C5.MGP"
meta[meta$fb_cluster==6,]$fb_celltype = "CAF.HLA"
meta[meta$fb_cluster==7,]$fb_celltype = "Fib.C7.IGKC"
meta[meta$fb_cluster==8,]$fb_celltype = "Fib.C8.SCN7A"
meta[meta$fb_cluster==9,]$fb_celltype = "Fib.C9.SCN7A"
meta[meta$fb_cluster==10,]$fb_celltype = "Fib.C10.PTGDS"
meta[meta$fb_cluster==11,]$fb_celltype = "Fib.C11.SLPI"
meta[meta$fb_cluster==12,]$fb_celltype = "Fib.C12.PLA2G2A"
meta[meta$fb_cluster==13,]$fb_celltype = "Fib.C13.PLA2G2A"
meta[meta$fb_cluster==14,]$fb_celltype = "CAF.ACTA2"
meta[meta$fb_cluster==15,]$fb_celltype = "Fib.C15.SCN7A"
meta[meta$fb_cluster==16,]$fb_celltype = "CAF.PDGFRB"
meta[meta$fb_cluster==17,]$fb_celltype = "Fib.C17.IGFBP2"
meta[meta$fb_cluster==18,]$fb_celltype = "CAF.CST1"
meta[meta$fb_cluster==19,]$fb_celltype = "CAF.PLA2G2A"
meta[meta$fb_cluster==20,]$fb_celltype = "Fib.C20.MYOC"
meta[meta$fb_cluster==21,]$fb_celltype = "Fib.C21.APOD"
meta[meta$fb_cluster==22,]$fb_celltype = "Fib.C22.PTGDS"
meta[meta$fb_cluster==23,]$fb_celltype = "Fib.C23.HSP"
meta[meta$fb_cluster==24,]$fb_celltype = "Fib.C24.PTGDS"
meta[meta$fb_cluster==25,]$fb_celltype = "CAF.MMP1"
meta[meta$fb_cluster==26,]$fb_celltype = "CAF.SPP1"
meta[meta$fb_cluster==27,]$fb_celltype = "CAF.TOP2A"
meta[meta$fb_cluster==28,]$fb_celltype = "Fib.C28.ACTA2"

#plot 
data@meta.data$celltype = meta$fb_celltype

pdf("figure/umap.cluster.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=FALSE, pt.size = 0.5)
dev.off()
pdf("figure/umap.celltype1.pdf",w=10,h=5)
DimPlot(data, reduction = "umap", group.by = "celltype", label=FALSE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("figure/umap.tissuetype.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "tissue", label=TRUE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#00BA38','Distant'='#619CFF'))
dev.off()
pdf("figure/umap.patient.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.5)
dev.off()

pdf("figure/PTGDS.pdf",w=7,h=5)
FeaturePlot(data, feature="PTGDS")
dev.off()
pdf("figure/SLPI.pdf",w=7,h=5)
FeaturePlot(data, feature="SLPI")
dev.off()
pdf("figure/PLA2G2A.pdf",w=7,h=5)
FeaturePlot(data, feature="PLA2G2A")
dev.off()
pdf("figure/SCN7A.pdf",w=7,h=5)
FeaturePlot(data, feature="SCN7A")
dev.off()
pdf("figure/MGP.pdf",w=7,h=5)
FeaturePlot(data, feature="MGP")
dev.off()
pdf("figure/MMP1.pdf",w=7,h=5)
FeaturePlot(data, feature="MMP1")
dev.off()
pdf("figure/CST1.pdf",w=7,h=5)
FeaturePlot(data, feature="CST1")
dev.off()
pdf("figure/COL1A1.pdf",w=7,h=5)
FeaturePlot(data, feature="COL1A1")
dev.off()
pdf("figure/ACTA2.pdf",w=7,h=5)
FeaturePlot(data, feature="ACTA2")
dev.off()
pdf("figure/SPP1.pdf",w=7,h=5)
FeaturePlot(data, feature="SPP1")
dev.off()
pdf("figure/VIM.pdf",w=7,h=5)
FeaturePlot(data, feature="VIM")
dev.off()
pdf("figure/SPARC.pdf",w=7,h=5)
FeaturePlot(data, feature="SPARC")
dev.off()
pdf("figure/FAP.pdf",w=7,h=5)
FeaturePlot(data, feature="FAP")
dev.off()
pdf("figure/PDGFRA.pdf",w=7,h=5)
FeaturePlot(data, feature="PDGFRA")
dev.off()
pdf("figure/PDGFRB.pdf",w=7,h=5)
FeaturePlot(data, feature="PDGFRB")
dev.off()


#gene identity
caf = meta[meta$FBcluster %in% c(7,11,12,13,21,22,24),]
data2 = data1[,rownames(caf)]

cd_gene = c("PTX3","MGP","SFRP1","MFAP5","SLPI","PLA2G2A","S100A10","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP7","TIMP1","TGM2","CXCL14","PDGFRA","COL1A1","COL3A1","ACTA2","POSTN","TAGLN","SPARC","CTHRC1","ACTB","PDGFRB","MMP1","MMP3","MMP13","MMP11","CXCL1","CXCL8","CST1","SDC1","S100A9","SRGN","C1QA","C1QB","HLA-DRA","HLA-DPA1","HLA-DPB1","FABP5","TOP2A","MKI67","CENPF","UBE2C")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))



fb = meta[meta$FBcluster %in% c(0,1,2,3,4,5,6,8,9,10,14,15,16,17,18,19,20,23),]
data3 = data1[,rownames(fb)]

cd_gene = c("SOD2","CXCL2","CCL2","PTGDS","PLA2G2A","RARRES1","MYOC","GSN","MFAP5","SLPI","S100A10","S100A4","IGFBP6","SCN7A","SPARCL1","PDGFRA","IGKC","IGFBP5","JUNB","SOCS3","SFRP4","HSPA1A","HSPA1B","HSP90AA1","IFI27","IFI6","CCL11","CCL13","IGFBP4","IGFBP7","VIM","CXCL14","TGFBI","PDGFRL","CYP1B1","DEPP1","APOE","GGT5","CFH","PTPRC","CXCR4","CCL5")
cd_gene = c("COL1A1","COL1A2","COL3A1","SPARC","VIM","PDGFRA","PDGFRB","ACTA2","FAP","TAGLN","CXCL8","POSTN","MMP1","MMP13","TOP2A","MKI67")

#vascular 

#dotplot
#1 CAF
caf = meta[meta$fb_cluster %in% c(6,14,16,18,19,25,26,27),]
data1 = data[,rownames(caf)]
cd_gene = c("MGP","PLA2G2A","S100A10","PTX3","SFRP1","SLPI","GSN","IGFBP6","ACTA2","COL1A1","COL3A1","POSTN","SPARC","SFRP4","MMP11","PDGFRB","CST1","IGFBP7","IGFBP2","IGFBP3","CXCL14","MMP1","MMP3","CXCL1","CXCL8","CXCL2","CCL2","SOD2","HLA-DRA","HLA-DPA1","HLA-DPB1","CXCR4","SRGN","C1QA","C1QB","APOE","SPP1","S100A9","FABP5","TOP2A","MKI67","CENPF","UBE2C")
data.usage = DotPlot(data1,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(19,14,16,18,25,6,26,27))
levels(data.usage$id) = c("CAF.PLG2A2G","CAF.ACTA2","CAF.PDGFRB","CAF.CST1","CAF.MMP1","CAF.HLA","CAF.SPP1","CAF.TOP2A")
pdf("figure/dotplot_CAF.pdf",h=4.5,w=17)
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

#2 normal
subset = meta[meta$fb_cluster %in% c(0,1,2,3,4,5,7,8,9,10,11,12,13,15,17,20,21,22,23,24,28),]
data1 = data[,rownames(subset)]
cd_gene = c("COL1A1","COL1A2","COL3A1","SPARC","VIM","PDGFRA","PDGFRB","ACTA2","FAP","TAGLN","CXCL8","POSTN","MMP1","MMP13","TOP2A","MKI67")
data.usage = DotPlot(data1,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(0,5,1,11,12,13,2,8,9,15,3,4,10,22,24,7,17,20,21,23,28))
levels(data.usage$id) = c("Fib.C0.MGP","Fib.C5.MGP","Fib.C1.SLPI","Fib.C11.SLPI","Fib.C12.PLA2G2A","Fib.C13.PLA2G2A","Fib.C2.PTX3","Fib.C8.SCN7A","Fib.C9.SCN7A","Fib.C15.SCN7A","Fib.C3.PTGDS","Fib.C4.PTGDS","Fib.C10.PTGDS","Fib.C22.PTGDS","Fib.C24.PTGDS","Fib.C7.IGKC","Fib.C17.IGFBP2","Fib.C20.MYOC","Fib.C21.APOD","Fib.C23.HSP","Fib.C28.ACTA2")
pdf("figure/dotplot_Normal.pdf",h=9,w=9)
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

#3 
subset = meta[meta$fb_cluster %in% c(0,1,2,3,4,5,7,8,9,10,11,12,13,15,17,20,21,22,23,24,28),]
data1 = data[,rownames(subset)]
cd_gene = c("MGP","DTP","CFH","GSN","SLPI","MFAP5","S100A4","S100A10","IGFBP6","PLA2G2A","RARRES1","PTX3","CXCL2","CCL2","SCN7A","SPARCL1","PDGFRA","PTGDS","APOE","APCDD1","A2M","CCL13","CCL11","IGFBP4","DEPP1","VIM","IGKC","IGFBP2","IGFBP7","MYOC","APOD","TGFBI","CYP1B1","HSPA1A","HSPA1B","HSP90AA1","ACTA2","TAGLN")
data.usage = DotPlot(data1,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(0,5,1,11,12,13,2,8,9,15,3,4,10,22,24,7,17,20,21,23,28))
levels(data.usage$id) = c("Fib.C0.MGP","Fib.C5.MGP","Fib.C1.SLPI","Fib.C11.SLPI","Fib.C12.PLA2G2A","Fib.C13.PLA2G2A","Fib.C2.PTX3","Fib.C8.SCN7A","Fib.C9.SCN7A","Fib.C15.SCN7A","Fib.C3.PTGDS","Fib.C4.PTGDS","Fib.C10.PTGDS","Fib.C22.PTGDS","Fib.C24.PTGDS","Fib.C7.IGKC","Fib.C17.IGFBP2","Fib.C20.MYOC","Fib.C21.APOD","Fib.C23.HSP","Fib.C28.ACTA2")
pdf("figure/dotplot_normalPheno.pdf",h=9,w=17)
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





#DotPlot
library(Seurat)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggplot2)
library(patchwork)
cd_gene = c("COL1A1","POSTN","TAGLN","SPARC","CXCL8","ACTA2","S100A4","VIM","FAP","PDGFRA","PDGFRB")
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

