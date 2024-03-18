##R 
##single cell data of ESCC
##step1: interpretation of the cell cluster
##Jiacheng Dai #daicy0424@gmail.com
##2023/02/02

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ROGUE)

#load data
data = readRDS("02data.23652.rds")
meta = read.csv("metadata.23652.csv",header=T,row.names=1)
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

#resolution = 1.3
#clustering
data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(object=data, resolution=1.3, reduction.type = "custom", k.param = 10, dims = 1:15, random.seed =888)
saveRDS(data1,"02data.res1.3.rds")

cluster.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers,file="02dataset_Markers.csv")

test1 = data1[[]]
meta$epi_cluster = test1$seurat_clusters
write.csv(meta,"metadata.res1.3.csv")

#celltype 
meta$epi_celltype = "NA"
meta[meta$epi_cluster == 0,]$epi_celltype = "Epi.N.C0"
meta[meta$epi_cluster == 1,]$epi_celltype = "Epi.T.C1"
meta[meta$epi_cluster == 2,]$epi_celltype = "Epi.T.C2"
meta[meta$epi_cluster == 3,]$epi_celltype = "Epi.T.C3"
meta[meta$epi_cluster == 4,]$epi_celltype = "Epi.T.C4"
meta[meta$epi_cluster == 5,]$epi_celltype = "Epi.T.C5"
meta[meta$epi_cluster == 6,]$epi_celltype = "Epi.T.C6"
meta[meta$epi_cluster == 7,]$epi_celltype = "Epi.T.C7"
meta[meta$epi_cluster == 8,]$epi_celltype = "Epi.T.C8"
meta[meta$epi_cluster == 9,]$epi_celltype = "Epi.N.C9"
meta[meta$epi_cluster == 10,]$epi_celltype = "Epi.N.C10"
meta[meta$epi_cluster == 11,]$epi_celltype = "Epi.T.C11"
meta[meta$epi_cluster == 12,]$epi_celltype = "Epi.T.C12"
meta[meta$epi_cluster == 13,]$epi_celltype = "Epi.T.C13"
meta[meta$epi_cluster == 14,]$epi_celltype = "Epi.T.C14"
meta[meta$epi_cluster == 15,]$epi_celltype = "Epi.T.C15"
meta[meta$epi_cluster == 16,]$epi_celltype = "Epi.T.C16"
meta[meta$epi_cluster == 17,]$epi_celltype = "Epi.N.C17"
meta[meta$epi_cluster == 18,]$epi_celltype = "Epi.T.C18"
meta[meta$epi_cluster == 19,]$epi_celltype = "Epi.N.C19"
meta[meta$epi_cluster == 20,]$epi_celltype = "Epi.T.C20"
meta[meta$epi_cluster == 21,]$epi_celltype = "Epi.T.C21"
meta[meta$epi_cluster == 22,]$epi_celltype = "Epi.N.C22"
meta[meta$epi_cluster == 23,]$epi_celltype = "Epi.N.C23"
meta[meta$epi_cluster == 24,]$epi_celltype = "Epi.N.C24"
meta[meta$epi_cluster == 25,]$epi_celltype = "Epi.T.C25"

data@meta.data$epi_celltype = meta$epi_celltype

#plot 
pdf("figure/umap_cluster.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("figure/umap_tissue.pdf",w=6,h=5)
DimPlot(data, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#00BA38','Distant'='#619CFF'))
dev.off()
pdf("figure/umap_patient.pdf",w=6,h=5)
DimPlot(data, reduction = "umap", group.by = "patient", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("figure/umap_celltype.pdf",w=8,h=5)
DimPlot(data, reduction = "umap", group.by = "epi_celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()


#Dotplot
library(ggplot2)
data = readRDS("02data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

meta1 = meta[meta$epi_cluster %in% c(1,2,3,4,5,6,7,8,11,12,13,14,15,16,18,20,21,25),]
data1 = data[,rownames(meta1)]
cd_gene = c("KRT5","SCPEP1","ITGA6","KRT4","ECM1","KRT15","EEF1A1","RPL3","RPS6","TFF3","MUC5B","MKI67","KRT7","KRT23")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/dist.pdf", h=5,w=5)
p
dev.off()




meta1 = meta[meta$epi_cluster %in% c(1,2,3,4,5,6,7,8,11,12,13,14,15,16,18,20,21,25),]
data1 = data[,rownames(meta1)]

cd_gene = c("AKR1C1","AKR1C2","AKR1C3","GSTM3","CES1","KRT14","LY6D","LGALS7B","CENPF","CENPW","TOP2A","HMGB1","PTTG1","MKI67","PDCD5")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p045_C1.pdf", h=5,w=8)
p
dev.off()

cd_gene = c("FGFBP2","CXCL14","FTL","CTSB","PTGR1","GPX2","KRT17","LY6E","IFI27")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p015_C7.pdf", h=4,w=8)
p
dev.off()

cd_gene = c("FABP5","S100A8","KRT6A","KRT6B","KRT14","KRT16","NDUFA4L2","NDRG1")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p045_C8.pdf", h=5,w=7)
p
dev.off()

cd_gene = c("COL1A1","COL1A2","COL3A1","MGP","DCN","VIM","EEF1A1","RPL3","RPS6")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p041_C14.pdf", h=4,w=8)
p
dev.off()


####################################################################################################################
cd_gene = c("COl1A1","COL1A2","COL3A1","MGP","DCN")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p049_C5.pdf", h=3,w=7)
p
dev.off()

cd_gene = c("SPRR1B","SPRR2A","S100A7","S100A9","SLPI","LCN2","CEACAM6","SPP1","VIM")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/dotplot/p049_C21.pdf", h=3,w=7)
p
dev.off()

meta1 = meta[meta$epi_cluster %in% c(10,13,31,17,29,19,28,23,14),]
data2 = data[,rownames(meta1)]

cd_gene = c("CEACAM6","PSCA","AGR2","LCN2","SAT1","GSTA1","KRT19","WFDC2","LY6E","KRT15","COL1A1","COL1A2","SPARC","TAGLN","CD74","C1QB","HLA-DQB1","SPRR2A","SPRR2D","S100A7","SLPI","MGP","COL3A1","VIM","SNAI2","ACKR1","PECAM1","KLK6","KLK7","PI3","SPP1","TPSB2","TPSAB1")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/Dotplot/p049_cplx.pdf", h=10,w=6)
p
dev.off()

meta1 = meta[meta$epi_cluster %in% c(4,5,7,15,22),]
data3 = data[,rownames(meta1)]

cd_gene = c("KRT13","KRT5","KRT15","AQP3","CRABP2","KRT4","KRT6C","CSTA","ANXA1","EEF1A1","RPL3","RPS6","CNFN","SPRR3","MAL","MGP","DCN","VIM","IGFBP4")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/Dotplot/p054_dist.pdf", h=10,w=6)
p
dev.off()


meta1 = meta[meta$epi_cluster %in% c(4,5,7,15,22),]
data3 = data[,rownames(meta1)]

cd_gene = c("KRT13","KRT5","KRT15","AQP3","CRABP2","KRT4","KRT6C","CSTA","ANXA1","EEF1A1","RPL3","RPS6","CNFN","SPRR3","MAL","MGP","DCN","VIM","IGFBP4")
p = DotPlot(data1, features = cd_gene) + coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
    scale_color_gradientn(values = seq(0,1,0.2), colors = c("#330066", "#336699", "#66CC66", "#FFCC33"))
pdf("figure/Dotplot/p054_dist.pdf", h=10,w=6)
p
dev.off()
#########################################################################################################################

##cell identity 
data = readRDS("02data.res1.3.rds")
expr = AverageExpression(data)
test1 = expr$RNA
test2 = expr$SCT

library(GSVA)
library(pheatmap)

genelist = read.csv("cellident1.csv",header=T,sep="\t")
colnames(genelist) = c("Celltype","Gene")
list = split(as.matrix(genelist)[,2], genelist[,1])

data = read.csv("expr.RNA.csv",header=T,row.names=1)
colnames(data) = paste0('C',0:25)
tumor = data[,c("C1","C2","C3","C4","C5","C6","C7","C8","C11","C12","C13","C14","C15","C16","C18","C20","C21","C25")]
dist = data[,c("C0","C9","C10","C17","C19","C22","C23","C24")]

data.gsva = gsva(as.matrix(tumor), list, method="ssgsea", kcdf = "Gaussian", abs.ranking=TRUE)
data.gsva1 = scale(data.gsva)
normalization = function(x) {
    return((x-min(x))/(max(x)-min(x)))
}
normdata.gsva = normalization(data.gsva1)
para = unique(c(seq(0,1,length=100)))

pheatmap(normdata.gsva,
         show_colnames = T,
         cluster_rows = F,cluster_cols = F,
         breaks=para)
