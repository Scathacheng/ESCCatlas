#R
#infercnv


library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

#cell annotation
setwd("03GeosAltas")

meta1 = read.csv("03tcell/metadata.res1.4.csv",header=T,row.names=1)
meta1.c4 = meta1[meta1$tcell_cluster == "4",]
res1 = meta1.c4[,c("tcell_cluster","tcell_celltype")]
colnames(res1) = c("sc_cluster","celltype")

meta = read.csv("02epithelial/metadata.res1.3.csv",header=T,row.names=1)
res2 = meta[,c("epi_cluster","epi_celltype")]
colnames(res2) = c("sc_cluster","celltype")

total = rbind(res1,res2)
write.csv(total, "02epithelial/infercnv/cellannotation.csv")

#data
matrix_dir="/home/public/project/singleCell/run/agg77/outs/filtered_feature_bc_matrix/"
barcode.path<-paste0(matrix_dir,"barcodes.tsv.gz")
features.path<-paste0(matrix_dir,"features.tsv.gz")
matrix.path<-paste0(matrix_dir,"matrix.mtx.gz")
data.matrix<-readMM(file=matrix.path)
feature.names=read.delim(features.path,header=F,stringsAsFactors=FALSE)
barcode.names=read.delim(barcode.path,header=F,stringsAsFactors=FALSE)
colnames(data.matrix)<-barcode.names$V1
rownames(data.matrix)<-feature.names$V2
data.matrix<-as.matrix(data.matrix)
data<-data.matrix[,rownames(total)]
saveRDS(data, "02epithelial/infercnv/infercnv_data.rds")


#analysis
library(Seurat)
library(infercnv)
data = readRDS("infercnv_data.rds")
cellanno<-read.table("cellannotation.txt",header=F,row.names=1,sep="\t")
cellanno<-as.matrix(cellanno)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=data,
					annotations_file=cellanno,
					delim="\t",
					gene_order_file="gencode_v21_gen_pos.txt",
					ref_group_names="TNK.C4")

infercnv_obj = infercnv::run(infercnv_obj,
			cutoff=0.1,
			out_dir="out",
			cluster_by_groups=TRUE,
			plot_steps=FALSE,
			denoise=TRUE,
			HMM=TRUE,
			no_prelim_plot=TRUE,
			output_format="pdf")

mcmc_obj = infercnv::inferCNVBayesNet(infercnv_obj = infercnv_obj,
			HMM_states = HMM_states,
			file_dir = "out",
			postMcmcMethod = "removeCNV",
			out_dir="out",
			resume_file_token="HMMi6.hmm_mode-samples",
			quietly = TRUE,
			CORES = 8,
			plotingProbs = FALSE,
			diagnostics = FALSE,
			HMM_type = "i6",
			k_obs_groups = 1,
			cluster_by_groups = FALSE,
			reassignCNVs = FALSE,
			no_plot = TRUE)

##CNV score
setwd("03GeosAltas/02epithelial/infercnv/")
object = read.table("out/infercnv.observations.txt",header=T)
object = as.matrix(object)
expr.scale = scale(t(object))

tmp1 = sweep(expr.scale, 2, apply(expr.scale, 2, min),"-")
tmp2 = apply(expr.scale, 2, max) - apply(expr.scale, 2, min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)
cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
rownames(cnv_score) = gsub("\\.","-",rownames(cnv_score))

meta = read.csv("../metadata.res1.3.csv",header=T,row.names=1)
meta$cnv_score = cnv_score[rownames(meta),]$cnv_score
write.csv(meta,"../metadata.res1.3.csv")


#plot
rm(list=ls())
data = readRDS("02data.res1.3.rds")
meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)


library(ggplot2)
library(ggpubr)
test = meta
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))
test$cnv_score = rescale(test$cnv_score, to=c(1,10))
pdf("figure/cnv/prop.pdf",h=4.5,w=4)
ggboxplot(test, x="tissue",y="cnv_score", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()

library(scales)
data@meta.data$cnv = rescale(meta$cnv_score, to=c(1,5))

pdf("figure/cnv/umap_cnv.pdf",w=5.5,h=5)
FeaturePlot(data, features="cnv", cols=c("#000080","#0000FF","#BEDFB5","#FFD700","#DC143C"))
dev.off()

#
