##myeloid cells

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ROGUE)

#load data
data = readRDS("06data.16636.rds")
meta = read.csv("metadata.16636.csv",header=T,row.names=1)
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


#res=0.9  all.rogue = 0.85
data = RunUMAP(data, reduction="harmony", dims = 1:30)
data = FindNeighbors(data, dims = 1:30) 
data = FindClusters(object=data, resolution=1.3)
saveRDS(data, "06data.res1.3.rds")

#Find cluster biomarkers
cluster.markers = FindAllMarkers(data,
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.5)
write.csv(cluster.markers, "06myeloid_markers.csv")


#
meta$myeloid_cluster = data@meta.data$seurat_clusters
meta$myeloid_celltype = "NA"
meta[meta$myeloid_cluster %in% c(0,1,5,19,21),]$myeloid_celltype = "M1-like Mac"
meta[meta$myeloid_cluster %in% c(4,6,7,9,10,12,16,18),]$myeloid_celltype = "M2-like Mac"
meta[meta$myeloid_cluster == 13,]$myeloid_celltype = "Monocyte"
meta[meta$myeloid_cluster %in% c(8,14,22),]$myeloid_celltype = "Neutrophil"
meta[meta$myeloid_cluster %in% c(2,20,3),]$myeloid_celltype = "cDC2"
meta[meta$myeloid_cluster == 15,]$myeloid_celltype = "cDC3"
meta[meta$myeloid_cluster %in% c(11,17),]$myeloid_celltype = "cDC1"
meta[meta$myeloid_cluster %in% c(25,26),]$myeloid_celltype = "proliferating myeloid"
meta[meta$myeloid_cluster %in% c(24,23),]$myeloid_celltype = "pDCs"

data@meta.data$celltype = meta$myeloid_celltype

pdf("figure/umap.cluster.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label=TRUE, pt.size = 0.5)
dev.off()
pdf("figure/umap.celltype.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "celltype", label=TRUE, pt.size = 0.5, repel=TRUE)
dev.off()
pdf("figure/umap.tissue.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "tissue", label=FALSE, pt.size = 0.5, cols=c('Tumor'='#F8766D','Adjacent'='#00BA38','Distant'='#619CFF'))
dev.off()
pdf("figure/umap.patient.pdf",w=7,h=5)
DimPlot(data, reduction = "umap", group.by = "patient", label=FALSE, pt.size = 0.5)
dev.off()

write.csv(meta, "metadata.res1.3.csv")

#cell identity
#all
meta1 = meta[meta$myeloid_cluster %in% c(0,6,13,8,3,15,11,24,26),]
data2 = data[,rownames(meta1)]
cd_gene = c("LILRA4","IL3RA","GZMB","CD1A","S100B","CLEC9A","FLT3","CD1C","FCER1A","CST3","AREG","ITGAX","IDO1","LAMP3","CCR7","S100A8","S100A9","FCGR3B","CSF3R","G0S2","CXCL8","CXCL1","C1QA","C1QB","C1QC","CD68","CD163","CCL4","CCL3","CD83","LYVE1","MKI67","TYMS","UBE2C")
data.usage = DotPlot(data2,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(24,11,3,15,13,8,6,0,26))
levels(data.usage$id) = c("pDCs","cDC1","cDC2","cDC3","Monocyte","Neutrophil","M2-like Macro","M1-like Macro","Proliferating myeloid")
pdf("figure/dotplot_all.pdf",h=5,w=16)
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


#DC
meta1 = meta[meta$myeloid_cluster %in% c(24,20,15,3,11),]
data2 = data[,rownames(meta1)]
cd_gene = c("LILRA4","GZMB","IL3RA","IRF4","IRF7","CXCR3","CXCR4","CLEC9A","S100B","CD1A","FLT3","IDO1","CADM1","CD1C","FCER1A","CD1E","AREG","NR4A3","ITGAX","IL1B","EREG","NLRP3","CCR7","LAMP3","CCL19","CCL22","FSCN1")
data.usage = DotPlot(data2,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(24,11,3,20,15))
levels(data.usage$id) = c("pDC.C24","cDC1.C11","cDC2.C3","cDC2.C20","cDC3.C15")
pdf("figure/dotplot_DC_HLA.pdf",h=5,w=12)
data.usage %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) +
geom_point() +
scale_size("% detected", range=c(0,10)) +
cowplot::theme_cowplot() +
theme(axis.line = element_blank()) +
theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=1)) +
ylab('') +
xlab('') +
theme(axis.ticks = element_blank()) +
scale_color_gradientn(colours=viridis::viridis(20), limits=c(-1,2.5))
dev.off()

cd_gene = c("BST1","F13A1","VCAN","S100A9","S100A8","RNASE2","FCN1","CD14","TMEM176B","PLBD1","MGST1","RAB3D","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5")

#Macrophage
meta1 = meta[meta$myeloid_cluster %in% c(21,1,5,9,6,12,16,0),]
data2 = data[,rownames(meta1)]
cd_gene = c("C1QA","C1QB","C1QC","MRC1","MERTK","CCL13","CCL18","CCL4","TNF","JAG1","CXCL10","CXCL11","CCL5","CD86","CXCL9","MMP9","MMP14","CD276","SPP1","CD44","CXCR4","VCAN","VEGFA","TNFAIP6")
data.usage = DotPlot(data2,features=cd_gene)$data
data.usage$id = factor(data.usage$id, levels=c(21,1,5,9,6,12,16,0))
levels(data.usage$id) = c("Macro.C21","Macro.C1","Macro.C5","Macro.C9","Macro.C6","Macro.C12","Macro.C16","Macro.C0")
pdf("figure/dotplot_Mac.pdf",h=5,w=12)
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

cd_gene = c("C1QA","C1QB","C1QC","LILRB5","LYVE1","CD163","CCL13","RNASE1","PLTP","FCGRT","CD86","APOE","APOC1","CTSD","FCGR3A","ISG15","SPP1","TNFAIP6","HLA-DPA1","HLA-DRA","CLEC10A","FCER1A","CHI3L1","CCL4","CCL3","CXCL2","TNF","VEGFA","NLRP3")


#计算phagocytosis score
library(ggpubr)
library(GSVA)

#c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
#c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")
#c("MRC1","CD163","MERTK","C1QA","C1QB","C1QC")
#c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
genelist = matrix(NA, nrow=81,ncol=2)
genelist[,1] = c(rep("M1",16),rep("M2",34),rep("Phagocytosis",6),rep("Angiogenesis",25))
genelist[,2] = c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7","IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4","MRC1","CD163","MERTK","C1QB","C1QA","C1QC","CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
genelist = as.data.frame(genelist)
colnames(genelist) = c("phenotype","gene")
genelist$phenotype = factor(genelist$phenotype, levels=c("M1","M2","Phagocytosis","Angiogenesis"))
list = split(as.matrix(genelist)[,2],genelist[,1])

data = readRDS("06data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

meta1 = meta[meta$myeloid_cluster %in% c(21,1,5,9,6,12,16,0),]
data1 = data[,rownames(meta1)]
data2 = AverageExpression(data1, group.by = "seurat_clusters")
expr = as.matrix(data2$SCT)
colnames(expr) = c("Macro.C0","Macro.C1","Macro.C5","Macro.C6","Macro.C9","Macro.C12","Macro.C16","Macro.C21")

data.gsva = gsva(as.matrix(expr), list, method="ssgsea", kcdf = "Gaussian", abs.ranking=TRUE)

test = t(data.gsva)
test = as.data.frame(test)
test$id = rownames(test)
test$id = factor(test$id, levels=c("Macro.C0","Macro.C16","Macro.C12","Macro.C6","Macro.C9","Macro.C5","Macro.C1","Macro.C21"))

normalization = function(x) {
    return((x-min(x))/(max(x)-min(x)))
}
scaledata = data.gsva
scaledata[1,] = normalization(data.gsva[1,])
scaledata[2,] = normalization(data.gsva[2,])
scaledata[3,] = normalization(data.gsva[3,])
scaledata[4,] = normalization(data.gsva[4,])
scaledata1 = t(scaledata)
scaledata1 = as.data.frame(scaledata1)
scaledata1$id = rownames(scaledata1)
pdf("figure/score.pdf",h=5,w=5)
ggscatter(scaledata1, x="Phagocytosis",y="Angiogenesis",color="id",label="id",repel=TRUE)
dev.off()


library(reshape2)
library(ggplot2)
library(ggpubr)
#myeloid cell-macrophage
setwd("03Geosaltas/06myeloidcell")
data = readRDS("06data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

mac = meta[meta$myeloid_celltype %in% c("M1-like Mac","M2-like Mac"),]
data1 = data[,rownames(mac)]
expr = data1@assays$SCT@counts

genelist = matrix(NA, nrow=81,ncol=2)
genelist[,1] = c(rep("M1",16),rep("M2",34),rep("Phagocytosis",6),rep("Angiogenesis",25))
genelist[,2] = c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7","IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4","MRC1","CD163","MERTK","C1QB","C1QA","C1QC","CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
genelist = as.data.frame(genelist)
colnames(genelist) = c("phenotype","gene")
genelist$phenotype = factor(genelist$phenotype, levels=c("M1","M2","Phagocytosis","Angiogenesis"))
list = split(as.matrix(genelist)[,2],genelist[,1])

data.gsva = gsva(as.matrix(expr),list,method = "ssgsea",kcdf = "Gaussian", abs.ranking=TRUE)
saveRDS(data.gsva, "MAC_GSVA.rds")

test = t(data.gsva)
test1 = as.data.frame(unlist(test))
test1$tissue = mac$tissue
m1 = mac[mac$myeloid_celltype == "M1-like Mac",]
m2 = mac[mac$myeloid_celltype == "M2-like Mac",]
m1data = test1[rownames(m1),]
m2data = test1[rownames(m2),]

test2 = melt(m2data)
test3 = test2[test2$variable=="Angiogenesis",]
test3$value = as.numeric(test3$value)
test3$tissue = factor(test3$tissue, levels = c("Tumor","Adjacent","Distant"))
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/angio_M2.pdf",h=5,w=4)
ggboxplot(test3, x="tissue",y="value", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()

#myeloid cell-DC
setwd("03Geosaltas/06myeloidcell")
data = readRDS("06data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

DC = meta[meta$myeloid_celltype %in% c("cDC1","cDC2","cDC3","pDCs"),]
data1 = data[,rownames(DC)]
expr = data1@assays$SCT@counts

genelist = matrix(NA,nrow=27,ncol=2)
genelist[,1] = c(rep("signature",3),rep("inflammatory",14),rep("antigenpresenting",10))
genelist[,2] = c("CD74","IL10RA","XCR1","BST1","CD163","F13A1","S100A9","S100A8","VCAN","RNASE2","FCN1","CD14","TMEM176B","PLBD1","MGST1","RAB3D","CD36","CD1C","HLA-A","HLA-B","HLA-C","HLA-DQA2","HLA-DQA1","HLA-DQB1","HLA-DPB1","HLA-DQB","FCGR2B")
genelist = as.data.frame(genelist)
colnames(genelist) = c("phenotype","gene")
genelist$phenotype = factor(genelist$phenotype, levels=c("signature","inflammatory","antigenpresenting"))
list = split(as.matrix(genelist)[,2],genelist[,1])

data.gsva = gsva(as.matrix(expr),list,method = "ssgsea",kcdf = "Gaussian", abs.ranking=TRUE)
saveRDS(data.gsva, "DC_GSVA.rds")

test = t(data.gsva)
test1 = as.data.frame(unlist(test))
test1$tissue = DC$tissue 
cdc1 = DC[DC$myeloid_celltype == "cDC1",]
cdc1data = test1[rownames(cdc1),]

test2 = melt(cdc1data)
test3 = test2[test2$variable=="antigenpresenting",]
test3$value = as.numeric(test3$value)
test3$tissue = factor(test3$tissue, levels = c("Tumor","Adjacent","Distant"))
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/ap_cdc1.pdf",h=5,w=4)
ggboxplot(test3, x="tissue",y="value", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()




#monocle
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(monocle)

setwd()
data = readRDS("06data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

#1 DC
meta1 = meta[meta$myeloid_cluster %in% c(11,17,2,20,3,15),]
meta1 = meta1[meta1$patient == "p045",]
data1 = data[,rownames(meta1)]
meta1$tissue = factor(meta1$tissue, levels=c("Tumor","Adjacent","Distant"))
#c(15,18,4,17,20)
#c("cDC1.CLEC9A+","cDC2.AREG-","cDC2.AREG+","cDC2.CD1A+","cDC3.LAMP3+")
data.matrix = as(as.matrix(data1@assays$RNA@counts),'sparseMatrix')
data.pd = new("AnnotatedDataFrame", data = meta1)
data.fd = data.frame(gene_short_name = row.names(data.matrix), row.names = row.names(data.matrix))
data.fd = new("AnnotatedDataFrame", data = data.fd)

hsmm <- newCellDataSet(data.matrix,
                       phenoData = data.pd,
                       featureData = data.fd,
                       expressionFamily = negbinomial.size(),
                       lowerDetectionLimit = 1) 

#QC: filtering low-quality cells
hsmm<-detectGenes(hsmm,min_expr=0.1)
print(head(fData(hsmm)))

expressed_genes<-row.names(subset(fData(hsmm),num_cells_expressed >= 5))
length(expressed_genes)

hsmm <- estimateSizeFactors(hsmm)
hsmm <- estimateDispersions(hsmm)

#trajectory step1: choose genes that define a cell's progress
diff_test_res<-differentialGeneTest(hsmm[expressed_genes,],fullModelFormulaStr="~percent.mt")
ordering_genes<-row.names(subset(diff_test_res,qval < 0.01))
hsmm<-setOrderingFilter(hsmm,ordering_genes)
#pdf("02gene used in step1.pdf")
#plot_ordering_genes(hsmm)
#dev.off()

#trajectory step2: reduce data dimensionality
hsmm<-reduceDimension(hsmm, max_components=2, method='DDRTree')

#trajectory step3: order cells along the trajectory
hsmm<-orderCells(hsmm)
pdf("monocle/monocleDC/02CellOrdersbyCluster.pdf")
plot_cell_trajectory(hsmm,color_by="myeloid_subtype")
dev.off()
pdf("monocle/monocleDC/02CellOrdersbycelltype.pdf")
plot_cell_trajectory(hsmm,color_by="myeloid_celltype")
dev.off()
pdf("monocle/monocleDC/02CellOrdersbyState.pdf")
plot_cell_trajectory(hsmm,color_by="State")
dev.off()
pdf("monocle/monocleDC/02CellOrdersbytissue.pdf")
plot_cell_trajectory(hsmm,color_by="tissue")
dev.off()
pdf("monocle/monocleDC/02CellOrdersbyPseudotime.pdf")
plot_cell_trajectory(hsmm,color_by="Pseudotime")
dev.off()

EC_state <- function(hsmm){
  if (length(unique(pData(hsmm)$State)) > 1){
    T3_counts <- table(pData(hsmm)$State, pData(hsmm)$myeloid_subtype)[,"cDC1.C11"]
    return(as.numeric(names(T3_counts)[which(T3_counts == max(T3_counts))]))
  }else {
    return (1)
  }
}
hsmm <- orderCells(hsmm, root_state = EC_state(hsmm))

#individual cluster's trajectory plot
pdf("monocle/monocleDC/02monocle2_cluster_indi.pdf",h=6,w=14)
plot_cell_trajectory(hsmm,color_by= "myeloid_subtype") +
    facet_wrap(~myeloid_subtype,nrow=2)
dev.off()
pdf("monocle/monocleDC/02monocle2_tissue_indi.pdf",h=6,w=14)
plot_cell_trajectory(hsmm,color_by= "tissue") +
    facet_wrap(~tissue,nrow=1)
dev.off()

#gene
my_genes = c("LAMP3","CCR7","CD1C","CLEC9A","XCR1","FCER1A","CD1A") 
hsmm_subset = hsmm[my_genes,]
plot_genes_in_pseudotime(hsmm_subset, color_by = "myeloid_celltype")
pdf("monocle/monocleDC/heatmap.pdf",h=5,w=4)
plot_pseudotime_heatmap(hsmm[my_genes,], num_clusters=2, cores=1, show_rownames=T)
dev.off()




#2 Macrophage
data = readRDS("06data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)

meta1 = meta[meta$myeloid_celltype %in% c("M1-like Mac","M2-like Mac"),]
meta1 = meta1[meta1$patient == "p055",]
data1 = data[,rownames(meta1)]
meta1$tissue = factor(meta1$tissue, levels=c("Tumor","Adjacent","Distant"))

data.matrix = as(as.matrix(data1@assays$RNA@counts),'sparseMatrix')
data.pd = new("AnnotatedDataFrame", data = meta1)
data.fd = data.frame(gene_short_name = row.names(data.matrix), row.names = row.names(data.matrix))
data.fd = new("AnnotatedDataFrame", data = data.fd)

hsmm <- newCellDataSet(data.matrix,
                       phenoData = data.pd,
                       featureData = data.fd,
                       expressionFamily = negbinomial.size(),
                       lowerDetectionLimit = 1) 

#QC: filtering low-quality cells
hsmm<-detectGenes(hsmm,min_expr=0.1)
print(head(fData(hsmm)))

expressed_genes<-row.names(subset(fData(hsmm),num_cells_expressed >= 5))
length(expressed_genes)

hsmm <- estimateSizeFactors(hsmm)
hsmm <- estimateDispersions(hsmm)

#trajectory step1: choose genes that define a cell's progress
diff_test_res<-differentialGeneTest(hsmm[expressed_genes,],fullModelFormulaStr="~percent.mt")
ordering_genes<-row.names(subset(diff_test_res,qval < 0.01))
hsmm<-setOrderingFilter(hsmm,ordering_genes)
#pdf("02gene used in step1.pdf")
#plot_ordering_genes(hsmm)
#dev.off()

#trajectory step2: reduce data dimensionality
hsmm<-reduceDimension(hsmm, max_components=2, method='DDRTree')

#trajectory step3: order cells along the trajectory
hsmm<-orderCells(hsmm)
pdf("monocle/monocleMac/02CellOrdersbyCluster.pdf")
plot_cell_trajectory(hsmm,color_by="myeloid_subtype")
dev.off()
pdf("monocle/monocleMac/02CellOrdersbyCelltype.pdf")
plot_cell_trajectory(hsmm,color_by="myeloid_celltype")
dev.off()
pdf("monocle/monocleMac/02CellOrdersbyState.pdf")
plot_cell_trajectory(hsmm,color_by="State")
dev.off()
pdf("monocle/monocleMac/02CellOrdersbytissue.pdf")
plot_cell_trajectory(hsmm,color_by="tissue")
dev.off()
pdf("monocle/monocleMac/02CellOrdersbyPseudotime.pdf")
plot_cell_trajectory(hsmm, color_by = "Pseudotime")
dev.off()

EC_state <- function(hsmm){
  if (length(unique(pData(hsmm)$State)) > 1){
    T3_counts <- table(pData(hsmm)$State, pData(hsmm)$myeloid_subtype)[,"Macro.C21"]
    return(as.numeric(names(T3_counts)[which(T3_counts == max(T3_counts))]))
  }else {
    return (1)
  }
}
hsmm <- orderCells(hsmm, root_state = EC_state(hsmm))

#individual cluster's trajectory plot
pdf("monocle/monocleMac/02monocle2_cluster_indi.pdf",h=6,w=10)
plot_cell_trajectory(hsmm,color_by= "myeloid_subtype") +
    facet_wrap(~myeloid_subtype,nrow=2)
dev.off()
pdf("monocle/monocleMac/02monocle2_tissue_indi.pdf",h=6,w=14)
plot_cell_trajectory(hsmm,color_by= "tissue") +
    facet_wrap(~tissue,nrow=1)
dev.off()

#gene
my_genes = c("C1QA","MRC1","MERTK","CCL13","CCL18","CD86","CD276","CCL4","CXCL10","CXCL9","SPP1","MMP14","MMP9") 
hsmm_subset = hsmm[my_genes,]
plot_genes_in_pseudotime(hsmm_subset, color_by = "myeloid_celltype")
pdf("monocle/monocleMac/heatmap.pdf",h=3,w=4)
plot_pseudotime_heatmap(hsmm[my_genes,], num_clusters=2, cores=1, show_rownames=T)
dev.off()



