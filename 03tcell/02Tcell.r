#plot T cell

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ROGUE)

#load data
data = readRDS("03data.47724.rds")
meta = read.csv("metadata.47724.csv",header=T,row.names=1)
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
            min.cell.n = 100)
    av.rogue <- c()
    for (j in 1:ncol(rogue.res)) {
    tmp.r <- rogue.res[,j]
    tmp.r <- tmp.r[!is.na(tmp.r)]
    av.rogue[j] <- mean(tmp.r)
    }
    all.rogue[i] = mean(na.omit(av.rogue))
}
all.rogue

#res = 2
#clustering 
data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(object=data, resolution=1.4, reduction = "custom", k.param = 10, dims = 1:15, random.seed = 1234)
saveRDS(data, "03data.res1.4.rds")

cluster.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers,file="03Tcell_Markers.csv")

test1 = data[[]]
meta$tcell_cluster = test1$seurat_clusters
write.csv(meta, "metadata.res1.4.csv")

data = readRDS()
pdf("figure/umap_cluster.pdf",w=8,h=8)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 1)
dev.off()

#DimPlot
data@meta.data$celltype = meta$Tcellsubtype
data@meta.data$specific = meta$tcell_celltype

pdf("figure/03tcell_cluster2.pdf",w=12,h=5)
DimPlot(data, reduction = "umap", group.by = "specific", label=FALSE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("figure/03tcell_tissue.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#00BA38','Distant'='#619CFF'))
dev.off()
pdf("figure/03tcell_patient.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.5)
dev.off()
pdf("figure/03tcell_celltype2.pdf",w=8,h=5)
DimPlot(data, reduction = "umap", group.by = "celltype", label=FALSE, pt.size = 0.5, repel=TRUE)
dev.off()


#celltype 
meta$Tcellsubtype="NA"
meta[meta$tcell_cluster == 0,]$Tcellsubtype = "naive CD8"
meta[meta$tcell_cluster == 1,]$Tcellsubtype = "naive CD8"
meta[meta$tcell_cluster == 2,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 3,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 4,]$Tcellsubtype = "cytotoxic T cell"
meta[meta$tcell_cluster == 5,]$Tcellsubtype = "DNAJB1 T cell"
meta[meta$tcell_cluster == 6,]$Tcellsubtype = "naive CD4"
meta[meta$tcell_cluster == 7,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 8,]$Tcellsubtype = "exhausted T cell"
meta[meta$tcell_cluster == 9,]$Tcellsubtype = "GZMK Tem"
meta[meta$tcell_cluster == 10,]$Tcellsubtype = "DNAJB1 T cell"
meta[meta$tcell_cluster == 11,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 12,]$Tcellsubtype = "exhausted T cell"
meta[meta$tcell_cluster == 13,]$Tcellsubtype = "IGFBP7 T cell"
meta[meta$tcell_cluster == 14,]$Tcellsubtype = "Treg"
meta[meta$tcell_cluster == 15,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 16,]$Tcellsubtype = "GZMK Tem"
meta[meta$tcell_cluster == 17,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 18,]$Tcellsubtype = "naive CD8"
meta[meta$tcell_cluster == 19,]$Tcellsubtype = "CXCR4 Tem"
meta[meta$tcell_cluster == 20,]$Tcellsubtype = "Treg"
meta[meta$tcell_cluster == 21,]$Tcellsubtype = "GZMK Tem"
meta[meta$tcell_cluster == 22,]$Tcellsubtype = "Treg"
meta[meta$tcell_cluster == 23,]$Tcellsubtype = "IGFBP7 T cell"
meta[meta$tcell_cluster == 24,]$Tcellsubtype = "cytotoxic T cell"
meta[meta$tcell_cluster == 25,]$Tcellsubtype = "naive CD8"
meta[meta$tcell_cluster == 26,]$Tcellsubtype = "naive T cell"
meta[meta$tcell_cluster == 27,]$Tcellsubtype = "cytotoxic T cell"
meta[meta$tcell_cluster == 28,]$Tcellsubtype = "cytotoxic T cell"
meta[meta$tcell_cluster == 29,]$Tcellsubtype = "DNAJB1 T cell"
meta[meta$tcell_cluster == 30,]$Tcellsubtype = "proliferating T cell"
meta[meta$tcell_cluster == 31,]$Tcellsubtype = "naive T cell"

meta$tcell_celltype = "NA"
meta[meta$tcell_cluster == 0,]$tcell_celltype = "Tn.CD8.C0"
meta[meta$tcell_cluster == 1,]$tcell_celltype = "Tn.CD8.C1"
meta[meta$tcell_cluster == 2,]$tcell_celltype = "Tn.C2"
meta[meta$tcell_cluster == 3,]$tcell_celltype = "Tn.C3"
meta[meta$tcell_cluster == 4,]$tcell_celltype = "Tcyt.C4"
meta[meta$tcell_cluster == 5,]$tcell_celltype = "Tcell.DNAJB1.C5"
meta[meta$tcell_cluster == 6,]$tcell_celltype = "Tn.CD4.C6"
meta[meta$tcell_cluster == 7,]$tcell_celltype = "Tn.C7"
meta[meta$tcell_cluster == 8,]$tcell_celltype = "Tex.C8"
meta[meta$tcell_cluster == 9,]$tcell_celltype = "Tem.GZMK.C9"
meta[meta$tcell_cluster == 10,]$tcell_celltype = "Tcell.DNAJB1.C10"
meta[meta$tcell_cluster == 11,]$tcell_celltype = "Tn.C11"
meta[meta$tcell_cluster == 12,]$tcell_celltype = "Tex.C12"
meta[meta$tcell_cluster == 13,]$tcell_celltype = "Tcell.IGFBP7.C13"
meta[meta$tcell_cluster == 14,]$tcell_celltype = "Treg.CCR6.C14"
meta[meta$tcell_cluster == 15,]$tcell_celltype = "Tn.C15"
meta[meta$tcell_cluster == 16,]$tcell_celltype = "Tem.GZMK.C16"
meta[meta$tcell_cluster == 17,]$tcell_celltype = "Tn.C17"
meta[meta$tcell_cluster == 18,]$tcell_celltype = "Tn.CD8.C18"
meta[meta$tcell_cluster == 19,]$tcell_celltype = "Tem.CXCR4.C19"
meta[meta$tcell_cluster == 20,]$tcell_celltype = "Treg.TNFRSF4.C20"
meta[meta$tcell_cluster == 21,]$tcell_celltype = "Tem.GZMK.C21"
meta[meta$tcell_cluster == 22,]$tcell_celltype = "Treg.CXCR3.C22"
meta[meta$tcell_cluster == 23,]$tcell_celltype = "Tcell.IGFBP7.C23"
meta[meta$tcell_cluster == 24,]$tcell_celltype = "Tcyt.C24"
meta[meta$tcell_cluster == 25,]$tcell_celltype = "Tn.CD8.C25"
meta[meta$tcell_cluster == 26,]$tcell_celltype = "Tn.C26"
meta[meta$tcell_cluster == 27,]$tcell_celltype = "Tcyt.C27"
meta[meta$tcell_cluster == 28,]$tcell_celltype = "Tem.GZMK.C28"
meta[meta$tcell_cluster == 29,]$tcell_celltype = "Tcell.DNAJB1.C29"
meta[meta$tcell_cluster == 30,]$tcell_celltype = "Tprof.CD4.C30"
meta[meta$tcell_cluster == 31,]$tcell_celltype = "Tn.C31"

#根据VlnPlot区分CD4和CD8细胞
pdf("figure/VlnPlot_CD3.pdf",h=5,w=14)
VlnPlot(data, features=c("CD3D","CD3E"), pt.size=0)
dev.off()

pdf("figure/VlnPlot_CD4.pdf",h=5,w=14)
VlnPlot(data, features=c("CD4","CD40LG"), pt.size=0)
dev.off()

pdf("figure/VlnPlot_CD8.pdf",h=5,w=14)
VlnPlot(data, features=c("CD8A","CD8B"), pt.size=0)
dev.off()

pdf("figure/FeaturePlot_CD4.pdf",w=7,h=5)
FeaturePlot(data, feature="CD4")
dev.off()

pdf("figure/FeaturePlot_CD8.pdf",w=7,h=5)
FeaturePlot(data, feature="CD8A")
dev.off()

pdf("figure/FeaturePlot_CD3.pdf",w=7,h=5)
FeaturePlot(data, feature="CD3D")
dev.off()

pdf("figure/FeaturePlot_GNLY.pdf",w=7,h=5)
FeaturePlot(data, feature="GNLY")
dev.off()

pdf("figure/FeaturePlot_GZMK.pdf",w=7,h=5)
FeaturePlot(data, feature="GZMK")
dev.off()

pdf("figure/FeaturePlot_PDCD1.pdf",w=7,h=5)
FeaturePlot(data, feature="PDCD1")
dev.off()



meta1 = meta[meta$tcell_cluster %in% c(2,0,16,19,24,12,10,13,6,14,30),]
data1 = data[,rownames(meta1)]
cd_gene = c("CD3D","CD3E","CD8A","CD8B","CCL5","GZMK","CXCR4","CCL4","GNLY","GZMB","PRF1","KLRD1","PDCD1","CXCL13","DNAJB1","HSPA1A","IGFBP7","CD4","CD40LG","EEF1A1","IL7R","CD27","CTLA4","TNFRSF4","FOXP3","MKI67","TOP2A",'UBE2C')
data.usage = DotPlot(data1, feature=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(2,0,16,19,24,12,10,13,6,14,30))
levels(data.usage$id) = c("naive T cell","CD8 T cell","Tem.GZMK","Tem.CXCR4","NK cell","exhausted T cell","Tcell.DNAJB1","Tcell.IGFBP7","CD4 T cell","Treg","proliferating CD4 T cell")
pdf("figure/dotplot.pdf",h=5,w=14)
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



#
library(ggplot2)

cd_gene = c("CD3D","CD4","FOXP3","IL2RA","IKZF2")
p = DotPlot(data2, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot_Treg.pdf", h=3,w=10)
p
dev.off()

cd_gene = c("CD3D","CD8A","IFNG")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot_cytotoxic.pdf", h=3,w=10)
p
dev.off()

cd_gene = c("CD3D","CD4","CD8A","PDCD1","HAVCR2","LAG3","CTLA4","TIGIT")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot_Tex.pdf", h=3,w=10)
p
dev.off()

cd_gene = c("CD3D","CD4","CD8A","TOP2A","MKI67","HMGB2","TUBA1B","UBE2C")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot_proliferate.pdf", h=3,w=10)
p
dev.off()

cd_gene = c("CD3D","CD8A","GZMA","GNLY","GZMB","GZMK","NKG7")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot_CD8effect.pdf", h=5,w=10)
p
dev.off()



cd_gene = c("CD3D","CD3E","CD3G","CD4","CD40LG","CD8A","CD8B","SELL","CCR7","NR4A1","RORC","TBX21","CCR6","CXCR3","IL17A","IFNG","CXCR5","PDCD1","FOXP3","CTLA4","SLC4A10","TRAV1-2","CX3CR1","GZMB","GZMH","GNLY","PRF1","CRTAM","GZMK","CCR9","ITGAE","ITGA1","SPRY1","SPRY2","TRDC","ITGAD","KIR2DL4","IKZF2","MKI67","TYROBP","NCAM1","FCGR3A","CD160","PCDH9","KIT","LST1")
cd_gene = c("TCF7","LEF1","CCR7","SELL","MAL")
cd_gene = c("IL7R","GPR183","ZFP36L2","CXCR4")
cd_gene = c("ZNF683","CD52","HOPX","ID2","CXCR6","XCL1","XCL2")
cd_gene = c("TBX21","ASCL2","CX3CR1","KLRG1")
cd_gene = c("KLRD1","TYROBP","KIR2DL3","KIR2DL1","KIR3DL1","KIR3DL2","CD160","EOMES","TXK","KLRC1","KIR2DL4")
cd_gene = c("GZMK","CXCR5","CCR4","CD28","CXCR3","GZMH","CD27","HLA-DRB1")
cd_gene = c("PDCD1","CXCL13","LAYN")
cd_gene = c("STAT1","IFIT1","ISG15","CCR1")
cd_gene = c("SLC4A10","KLRB1","TMIGD2","RORA","RORC","ZBTB16","IL26","IL17A","IL23R")
cd_gene = c("NME1","NME2","MND1","SPC24","MYB")
cd_gene = c("RTKN2","IL2RA","S1PR1","TNFRSF9","CTLA4","LAYN","IFIT1","IRF7")

#分CD8 和 CD4的细胞分别画图
meta.cd4 = meta1[meta1$seurat_clusters %in% c(1,3,5,6,18,14,16),]
data.cd4 = data1[,rownames(meta.cd4)]

meta.cd8 = meta1[meta1$seurat_clusters %in% c(0,2,4,5,7,8,9,10,11,12,13,14,15,16,17,19),]
data.cd8 = data1[,rownames(meta.cd8)]

cd_gene = c("CD3D","CD3E","CD3G","CD4","CD40LG","SELL","FOXP3","CTLA4","CCR6","IKZF2","CCR7","IL7R","TNFRSF9","LAYN","MKI67")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot.CD4.Treg.pdf",h=5,w=6)
p
dev.off()

cd_gene = c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMB","GZMH","GNLY","PRF1","CRTAM","GZMK","TBX21","CX3CR1","ZNF683","CD52","CXCR4","CXCR6")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot.CD8.Tem.pdf",h=4,w=6)
p
dev.off()

cd_gene = c("CD3D","CD3E","CD3G","CD8A","CD8B","PDCD1","CTLA4","CXCR6","CXCL13","LAYN")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot.CD8.Texhaust.pdf",h=3.5,w=6)
p
dev.off()

cd_gene = c("CD3D","CD3E","CD3G","CD4","CD40LG","CD8A","CD8B","IL7R","ZFP36L2","IFNG","GZMB","GNLY","TRDC","TYROBP","HOPX","ID2","XCL1","XCL2","KLRD1","KLRC1")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/DotPlot.CD8.NK.pdf",h=6,w=6)
p
dev.off()




###计算score
data = readRDS("03data.res1.4.rds")
data1 = AverageExpression(object=data, group.by="seurat_clusters")
expr = as.matrix(data1$SCT)
colnames(expr) =  paste0('C',0:31)

#c("GZMB","PRF1")
#c("ZWINT","E2F1","FEN1","FOXM1","H2AFZ","HMGB2","MCM2","MCM3","MCM4","MCM5","MCM6","MKI67","MYBL2","PCNA","PLK1","CCND1","AURKA","BUB1","TOP2A","TYMS","DEK","CCNB1","CCNE1")

genelist = matrix(NA, nrow=25,ncol=2)
genelist[,1] = c(rep("cytotoxic",2),rep("proliferate",23))
genelist[,2] = c("GZMB","PRF1","ZWINT","E2F1","FEN1","FOXM1","H2AFZ","HMGB2","MCM2","MCM3","MCM4","MCM5","MCM6","MKI67","MYBL2","PCNA","PLK1","CCND1","AURKA","BUB1","TOP2A","TYMS","DEK","CCNB1","CCNE1")
genelist = as.data.frame(genelist)
colnames(genelist) = c("phenotype","gene")
list = split(as.matrix(genelist)[,2],genelist[,1])

data.gsva = gsva(as.matrix(expr), list, method="ssgsea", kcdf = "Gaussian", abs.ranking=TRUE)

cy = expr[c("GZMB","PRF1"),]
for (i in 1:32){
    data.gsva[1,i] = sqrt(cy[1,i]*cy[2,i])
}
data.gsva[1,] = log10(data.gsva[1,]+1)

library(ggpubr)
test = t(data.gsva)
test = as.data.frame(test)
test$id = rownames(test)
test$id = factor(test$id, levels=c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31"))
pdf("figure/score.pdf",h=7,w=5)
ggscatter(test, x="cytotoxic",y="proliferate",color="id", label="id",repel=TRUE)
dev.off()


###marker gene
#1 CD4+
meta1 = meta[meta$tcell_cluster %in% c(2,6,14,22,20,30),]
data2 = data[,rownames(meta1)]

cd_gene = c("CD3D","CD3E","CD4","CD40LG","EEF1A1","IL7R","IKZF2","CD27","CTLA4","FOXP3","TIGIT","BATF","IL32","CXCR3","CCR6","TNFRSF4","TNFRSF9","TNFRSF18","TOP2A","MKI67","UBE2C")
p = DotPlot(data2, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))

pdf("figure/CD4.pdf",h=4,w=10)
data.usage = DotPlot(data2, feature=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(2,6,22,14,20,30))
levels(data.usage$id) = c("naive T cell","Tn.CD4.C6","Treg.C22","Treg.C14","Treg.C20","Tprof.C30")
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



#2 effect T cell
meta1 = meta[meta$tcell_cluster %in% c(1,9,16,19,28),]
data2 = data[,rownames(meta1)]

cd_gene = c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMK","CCL5","CCL4","CXCR4","IFNG","NKG7","GNLY","GZMB","PRF1","TBX21","ASCL2","ZNF683","CD52","GZMH","CXCL13","PDCD1","CTLA4","TIGIT")
cd_gene = c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMK","CCL5","GZMH","NKG7","GZMA","CST7","CCL4","CXCR4","CCL4L2","CCL3L1","GNLY","GZMB","PRF1")

p = DotPlot(data2, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))

pdf("figure/Tem.pdf",h=4,w=10)
data.usage = DotPlot(data2, feature=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(1,9,16,28,19))
levels(data.usage$id) = c("Tn.CD8.C0","Tem.GZMK.C9","Tem.GZMK.C16","Tem.GZMK.C28","Tem.CXCR4.C19")
data.usage %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5)) +
theme(axis.title.x = element_blank()) +
ylab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(0,2))
dev.off()

#3 exhausted 
meta1 = meta[meta$tcell_cluster %in% c(1,8,12,4,24,27,28),]
data2 = data[,rownames(meta1)]

cd_gene = c("CD3D","CD3E","CD3G","CD8A","CD8B","GNLY","GZMB","PRF1","KLRD1","KLRC1","LAG3","HAVCR2","TIGIT","PDCD1","CTLA4","CXCL13")

p = DotPlot(data2, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))

pdf("figure/Tex.pdf",h=4,w=10)
data.usage = DotPlot(data2, feature=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(1,24,27,4,28,8,12))
levels(data.usage$id) = c("Tn.CD8.C0","Tcyt.C24","Tcyt.C27","Tcyt.C4","Tem.GZMK.C28","Tex.C8","Tex.C12")
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

#NK
meta1 = meta[meta$Tcell_cluster %in% c(3,5,7,8,10,11,16,17,19),]
data2 = data[,rownames(meta1)]

cd_gene = c("CD3D","CD3E","CD3G","CD4","CD40LG","CD8A","CD8B","IL7R","ZFP36L2","KIR2DL3","KIR2DL1","KIR3DL1","KIR3DL2","GZMB","GNLY","TRDC","TYROBP","HOPX","ID2","XCL1","XCL2","KLRD1","KLRC1")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))


#other
meta1 = meta[meta$Tcell_cluster %in% c(3,5,7,8,10,16,17,19),]
data2 = data[,rownames(meta1)]

cd_gene = c("CD3D","CD4","CD8A","DNAJB1","TXNIP","CCR7","SELL","RPS3","RPL6","SPARC","IGFBP7","IGKC","XIST","SYNE2","IRF7","STAT1","IFIT1","ISG15","CCR1","MS4A1","CD79A","CD74","HLA-DRA","HSPA1A")
p = DotPlot(data, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))

##T cell-different region diff
library(reshape2)
library(ggplot2)
library(ggpubr)
setwd("03Geosaltas/03tcell")
data=readRDS("03data.res1.4.rds")
meta=read.csv("metadata.res1.4.csv",header=T,row.names=1)
expr = data@assays$SCT@counts


genelist = matrix(NA, nrow=36,ncol=2)
genelist[,1] = c(rep("cytotoxic",2),rep("cytotoxic2",2),rep("proliferate",23),rep("treg",9))
genelist[,2] = c("GZMB","PRF1","NKG7","GNLY","ZWINT","E2F1","FEN1","FOXM1","H2AFZ","HMGB2","MCM2","MCM3","MCM4","MCM5","MCM6","MKI67","MYBL2","PCNA","PLK1","CCND1","AURKA","BUB1","TOP2A","TYMS","DEK","CCNB1","CCNE1","CD27","IKZF2","BATF","IL32","CXCR3","CCR6","TNFRSF4","TNFRSF9","TNFRSF18")
genelist = as.data.frame(genelist)
colnames(genelist) = c("phenotype","gene")
list = split(as.matrix(genelist)[,2],genelist[,1])

data.gsva = gsva(as.matrix(expr),list,method = "ssgsea",kcdf = "Gaussian", abs.ranking=TRUE)
saveRDS(data.gsva, "T_GSVA.rds")

test = t(data.gsva)
test1 = as.data.frame(unlist(test))
test1$tissue = meta$tissue

subtype = meta[meta$Tcellsubtype %in% c("Tem.GZMK","Tem.CXCR4"),]
test2 = test1[rownames(subtype),]
test2 = melt(test2)

#cytotoxic
#cytotoxic2
#proliferate
#treg

test3 = test2[test2$variable=="cytotoxic2",]
test3$value = as.numeric(test3$value)
test3$tissue = factor(test3$tissue, levels = c("Tumor","Adjacent","Distant"))
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/cyto2_Tem.pdf",h=5,w=4)
ggboxplot(test3, x="tissue",y="value", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()



