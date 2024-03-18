#infercnv

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

#infercnv cellannotation
#example
#p045
rm(list=ls())
meta1 = read.csv("../new0418/metadata.211.csv",header=T,row.names=1)
meta1.c3 = meta1[meta1$Patient_ID == "p045" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2 = read.csv("02Metadata_res1.1.csv",header=T,row.names=1)
meta2.p045 = meta2[meta2$Patient_ID == "p045",]
meta2.p045 = meta2.p045[meta2.p045$epi_cluster %in% c(0,1,2,3,6,11,16,17,20,22,26),] #which lager than 100
meta2.res = meta2.p045[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p045/cellannotation.csv")
#need edit
#remove the column of Patient_ID
#remove the colname

#p026
meta1.c3 = meta1[meta1$Patient_ID == "p026" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p026 = meta2[meta2$Patient_ID == "p026",]
meta2.p026 = meta2.p026[meta2.p026$epi_cluster %in% c(4,6,16),] #which lager than 100
meta2.res = meta2.p026[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p026/cellannotation.csv")


#p036
meta1.c3 = meta1[meta1$Patient_ID == "p036" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p036 = meta2[meta2$Patient_ID == "p036",]
meta2.p036 = meta2.p036[meta2.p036$epi_cluster %in% c(4,5,6,7,15,16,17,19,21,22,26),] #which lager than 100
meta2.res = meta2.p036[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p036/cellannotation.csv")


#p037
meta1.c3 = meta1[meta1$Patient_ID == "p037" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p037 = meta2[meta2$Patient_ID == "p037",]
meta2.p037 = meta2.p037[meta2.p037$epi_cluster %in% c(4,7,21,25),] #which lager than 100
meta2.res = meta2.p037[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p037/cellannotation.csv")


#p041
meta1.c3 = meta1[meta1$Patient_ID == "p041" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p041 = meta2[meta2$Patient_ID == "p041",]
meta2.p041 = meta2.p041[meta2.p041$epi_cluster %in% c(9,15,21,26),] #which lager than 100
meta2.res = meta2.p041[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p041/cellannotation.csv")


#p047
meta1.c3 = meta1[meta1$Patient_ID == "p047" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p047 = meta2[meta2$Patient_ID == "p047",]
meta2.p047 = meta2.p047[meta2.p047$epi_cluster %in% c(18,26,27),] #which lager than 100
meta2.res = meta2.p047[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p047/cellannotation.csv")


#p048
meta1.c3 = meta1[meta1$Patient_ID == "p048" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p048 = meta2[meta2$Patient_ID == "p048",]
meta2.p048 = meta2.p048[meta2.p048$epi_cluster %in% c(13,24),] #which lager than 100
meta2.res = meta2.p048[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p048/cellannotation.csv")


#p049
meta1.c3 = meta1[meta1$Patient_ID == "p049" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p049 = meta2[meta2$Patient_ID == "p049",]
meta2.p049 = meta2.p049[meta2.p049$epi_cluster %in% c(10,13,14,17,19,23,28,29,30,31),] #which lager than 100
meta2.res = meta2.p049[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p049/cellannotation.csv")


#p050
meta1.c3 = meta1[meta1$Patient_ID == "p050" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p050 = meta2[meta2$Patient_ID == "p050",]
meta2.p050 = meta2.p050[meta2.p050$epi_cluster %in% c(8,7,13),] #which lager than 100
meta2.res = meta2.p050[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p050/cellannotation.csv")



#p054
meta1.c3 = meta1[meta1$Patient_ID == "p054" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p054 = meta2[meta2$Patient_ID == "p054",]
meta2.p054 = meta2.p054[meta2.p054$epi_cluster %in% c(4,5,7,13,15,16,22),] #which lager than 100
meta2.res = meta2.p054[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p054/cellannotation.csv")


#p055
meta1.c3 = meta1[meta1$Patient_ID == "p055" & meta1$cluster == "3",]
meta1.res = meta1.c3[,c("Patient_ID","cluster")]
colnames(meta1.res) = c("Patient_ID","seurat_clusters")
meta1.res[,2] = "NK_C3"
meta2.p055 = meta2[meta2$Patient_ID == "p055",]
meta2.p055 = meta2.p055[meta2.p055$epi_cluster %in% c(4,6,7,12,22),] #which lager than 100
meta2.res = meta2.p055[,c("Patient_ID","epi_cluster")]
colnames(meta2.res) = c("Patient_ID","seurat_clusters")
meta.res = rbind(meta1.res, meta2.res)
write.csv(meta.res, "./infercnv/p055/cellannotation.csv")





#on server
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
data<-data.matrix[,rownames(meta.res)]
saveRDS(data, "./infercnv/p045/infercnv.rds")

#on PC
library(Seurat)
library(infercnv)
data = readRDS("./p055/infercnv.rds")
cellanno<-read.table("./p055/cellannotation.txt",header=F,row.names=1,sep="\t")
cellanno<-as.matrix(cellanno)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=data,
					annotations_file=cellanno,
					delim="\t",
					gene_order_file="gencode_v21_gen_pos.txt",
					ref_group_names="TNK.C4")

infercnv_obj = infercnv::run(infercnv_obj,
			cutoff=0.1,
			out_dir="p055",
			cluster_by_groups=TRUE,
			plot_steps=FALSE,
			denoise=TRUE,
			HMM=TRUE,
			no_prelim_plot=TRUE,
			output_format="pdf")

mcmc_obj = infercnv::inferCNVBayesNet(infercnv_obj = infercnv_obj,
			HMM_states = HMM_states,
			file_dir = "p055",
			postMcmcMethod = "removeCNV",
			out_dir="p055",
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



object = readRDS("test3/p055/17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj")
cnv_value = object@expr.data
cnv_value = (cnv_value-3)/2
cnv_level = colMeans(abs(cnv_value))
table(cnv_level) #选择一个阈值讨论
cellanno<-read.csv("test3/p055/cellannotation.csv",header=T,row.names=1)
table(cellanno)
cellanno$cnv = cnv_level
p055 = cellanno


result = rbind(p026, p036, p037, p041, p045, p047, p048, p049, p050, p054, p055)
res1 = table(x=result$seurat_clusters, y=result$cnv)
res1 = as.data.frame(res1)
res1 = res1[res1$Freq > 0, ]
colnames(res1)=c("SeuratCluster","CNV","Freq")
write.csv(res1,"CNV_res1.csv")


meta = read.csv("02Metadata_res1.1.csv",header=T,row.names=1)
result = read.csv("total_cnv.csv",header=T,row.names=1)
result$tissue = meta[match(rownames(result),rownames(meta)),]$tissue

library(ggplot2)
library(ggpubr)
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))
ggboxplot(meta, x="tissue",y="cnv_score", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)






