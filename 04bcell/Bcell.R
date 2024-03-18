# select B cell and plasma cell

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ROGUE)

#load data
data = readRDS("04data.33202.rds")
meta = read.csv("metadata.33202.csv",header=T,row.names=1)
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

data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(object=data, resolution=1.1)
saveRDS(data, "04data.res1.1.rds")

cluster.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers,file="04Bcell_Markers.csv")

test1 = data[[]]
meta$bcell_cluster = test1$seurat_clusters
write.csv(meta, "metadata.res1.1.csv")

#identity
meta$bcell_celltype = "NA"
meta[meta$bcell_cluster == 0,]$bcell_celltype = "Mem B1"
meta[meta$bcell_cluster == 1,]$bcell_celltype = "Mem B2"
meta[meta$bcell_cluster == 2,]$bcell_celltype = "Mem B3"
meta[meta$bcell_cluster == 3,]$bcell_celltype = "Stress B1"
meta[meta$bcell_cluster == 4,]$bcell_celltype = "Mem B4"
meta[meta$bcell_cluster == 5,]$bcell_celltype = "IgG3 PC3"
meta[meta$bcell_cluster == 6,]$bcell_celltype = "Mem B5"
meta[meta$bcell_cluster == 7,]$bcell_celltype = "Mem B6"
meta[meta$bcell_cluster == 8,]$bcell_celltype = "Mem B7"
meta[meta$bcell_cluster == 9,]$bcell_celltype = "IgA PC1"
meta[meta$bcell_cluster == 10,]$bcell_celltype = "IgG3 PC9"
meta[meta$bcell_cluster == 11,]$bcell_celltype = "IgA PC2"
meta[meta$bcell_cluster == 12,]$bcell_celltype = "IgA PC3"
meta[meta$bcell_cluster == 13,]$bcell_celltype = "IgA PC4"
meta[meta$bcell_cluster == 14,]$bcell_celltype = "IgG3 PC8"
meta[meta$bcell_cluster == 15,]$bcell_celltype = "IgG3 PC1"
meta[meta$bcell_cluster == 16,]$bcell_celltype = "IgA PC5"
meta[meta$bcell_cluster == 17,]$bcell_celltype = "IgG4 PC2"
meta[meta$bcell_cluster == 18,]$bcell_celltype = "IgG4 PC1"
meta[meta$bcell_cluster == 19,]$bcell_celltype = "IgG2 PC1"
meta[meta$bcell_cluster == 20,]$bcell_celltype = "IgG3 PC4"
meta[meta$bcell_cluster == 21,]$bcell_celltype = "IgG3 PC2"
meta[meta$bcell_cluster == 22,]$bcell_celltype = "IgG3 PC5"
meta[meta$bcell_cluster == 23,]$bcell_celltype = "IgG2 PC2"
meta[meta$bcell_cluster == 24,]$bcell_celltype = "IgA PC6"
meta[meta$bcell_cluster == 25,]$bcell_celltype = "IgG2 PC3"
meta[meta$bcell_cluster == 26,]$bcell_celltype = "IgG3 PC6"
meta[meta$bcell_cluster == 27,]$bcell_celltype = "Naive B cell"
meta[meta$bcell_cluster == 28,]$bcell_celltype = "IgG2 PC4"
meta[meta$bcell_cluster == 29,]$bcell_celltype = "IgM PC"
meta[meta$bcell_cluster == 30,]$bcell_celltype = "IgA PC7"
meta[meta$bcell_cluster == 31,]$bcell_celltype = "IgA PC8"
meta[meta$bcell_cluster == 32,]$bcell_celltype = "Stress B2"
meta[meta$bcell_cluster == 33,]$bcell_celltype = "Mem B8"
meta[meta$bcell_cluster == 34,]$bcell_celltype = "IgG3 PC7"


data@meta.data$celltype = meta$bcell_celltype
pdf("figure/umap_cluster.pdf",w=8,h=6)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 1)
dev.off()
pdf("figure/umap_tissue.pdf",w=8,h=6)
DimPlot(data, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.1, cols=c('Tumor'='#F8766D','Adjacent'='#00BA38','Distant'='#619CFF'))
dev.off()
pdf("figure/umap_patient.pdf",w=8,h=6)
DimPlot(data, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 1)
dev.off()
pdf("figure/umap_celltype.pdf",w=10,h=6)
DimPlot(data, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 1, repel = TRUE)
dev.off()

cd_gene = c("MS4A1","CD19","CXCR4","CD37","CD79B","CD24","BANK1","SELL","HSPA1A","HSPA6","IL4R","TCL1A","BACH2","MKI67","TUBB","STMN1","TYMS","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","IGKC","IGLC2","IGLC3","CD27","CD38","JCHAIN","MZB1","SDC1","TNFRSF17")
data.usage = DotPlot(data,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(27,0,1,2,4,6,7,8,33,3,32,9,11,12,13,16,24,30,31,19,23,25,28,15,21,5,20,22,26,34,14,10,18,17,29))
levels(data.usage$id) = c("Naive B","Mem B1","Mem B2","Mem B3","Mem B4","Mem B5","Mem B6","Mem B7","Mem B8","Stress B1","Stress B2","IgA PC1","IgA PC2","IgA PC3","IgA PC4","IgA PC5","IgA PC6","IgA PC7","IgA PC8","IgG2 PC1","IgG2 PC2","IgG2 PC3","IgG2 PC4","IgG3 PC1","IgG3 PC2","IgG3 PC3","IgG3 PC4","IgG3 PC5","IgG3 PC6","IgG3 PC7","IgG3 PC8","IgG3 PC9","IgG4 PC1","IgG4 PC2","IgM PC")
pdf("figure/dotplot.pdf",h=12,w=15)
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




