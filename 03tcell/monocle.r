library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(monocle)

setwd()
data = readRDS("06data.res1.4.rds")
meta = read.csv("metadata.res1.4.csv",header=T,row.names=1)

#1 DC
meta1 = meta[meta$tcell_cluster %in% c(8,12,24,4,0,9,21,16,19,27,28,1),]
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
diff_test_res<-differentialGeneTest(hsmm[expressed_genes,],fullModelFormulaStr="~tcell_cluster")
ordering_genes<-row.names(subset(diff_test_res,qval < 0.01))
hsmm<-setOrderingFilter(hsmm,ordering_genes)
#pdf("02gene used in step1.pdf")
#plot_ordering_genes(hsmm)
#dev.off()

#trajectory step2: reduce data dimensionality
hsmm<-reduceDimension(hsmm, max_components=2, method='DDRTree')

#trajectory step3: order cells along the trajectory
hsmm<-orderCells(hsmm)
pdf("monocle/02CellOrdersbyCluster.pdf")
plot_cell_trajectory(hsmm,color_by="tcell_celltype")
dev.off()
pdf("monocle/02CellOrdersbyState.pdf")
plot_cell_trajectory(hsmm,color_by="State")
dev.off()
pdf("monocle/02CellOrdersbytissue.pdf")
plot_cell_trajectory(hsmm,color_by="tissue")
dev.off()

EC_state <- function(hsmm){
  if (length(unique(pData(hsmm)$State)) > 1){
    T3_counts <- table(pData(hsmm)$State, pData(hsmm)$tcell_celltype)[,"Tn.CD8.C0"]
    return(as.numeric(names(T3_counts)[which(T3_counts == max(T3_counts))]))
  }else {
    return (1)
  }
}
hsmm <- orderCells(hsmm, root_state = EC_state(hsmm))
pdf("monocle/02monocle2_pseudotime.pdf")
plot_cell_trajectory(hsmm, color_by = "Pseudotime")
dev.off()
pdf("monocle/02monocle2_cluster.pdf")
plot_cell_trajectory(hsmm, color_by = "tcell_celltype", pt.size = .1)
dev.off()
pdf("monocle/02monocle2_tissue.pdf")
plot_cell_trajectory(hsmm, color_by = "tissue")
dev.off()

#individual cluster's trajectory plot
pdf("monocle/02monocle2_cluster_indi.pdf",h=6,w=14)
plot_cell_trajectory(hsmm,color_by= "tcell_celltype") +
    facet_wrap(~tcell_celltype,nrow=2)
dev.off()
pdf("monocle/02monocle2_tissue_indi.pdf",h=6,w=14)
plot_cell_trajectory(hsmm,color_by= "tissue") +
    facet_wrap(~tissue,nrow=1)
dev.off()


