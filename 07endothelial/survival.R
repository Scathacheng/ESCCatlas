#survival analysis in ESCC
library(survival)
library(survminer)

#load data
setwd("Downloads/04survival/")
meta = read.table("lusc/TCGA-LUSC.survival.tsv",header=T,row.names=1,sep="\t",quote="")
#meta1 = read.table("esca/TCGA-ESCA.GDC_phenotype.tsv",header=T,row.names=1,sep="\t",quote="")
#meta$disease = meta1[rownames(meta),]$disease_type
#meta = meta[meta$disease=="Squamous Cell Neoplasms",]
rownames(meta) = gsub("-",".",rownames(meta))

matrix = read.table("lusc/TCGA-LUSC.htseq_counts.tsv.gz",header=T,row.names=1)
list = read.csv("gencode.v22.annotation.gene.probeMap.csv",header=T,sep="\t",row.names=1)
geneoverlap = intersect(rownames(list),rownames(matrix))
sampleoverlap = intersect(rownames(meta),colnames(matrix))

matrix = matrix[geneoverlap,sampleoverlap]
list = list[geneoverlap,]
meta.survival = meta[sampleoverlap,]
#esca 82
#luad 572
#lusc 542
#change gene id
matrix$gene = list$gene
matrix1 = matrix[!duplicated(matrix$gene),]
rownames(matrix1) = matrix1[,543]
matrix1 = matrix1[,-543]
expr = matrix[,1:542]
s = unique(matrix$gene[duplicated(matrix$gene)])
matrix2 = sapply(s, function(i){
	colMeans(expr[matrix$gene==i,])
	})
matrix2 = t(matrix2)
test = matrix1[setdiff(rownames(matrix1),s),]
matrix3 = rbind(matrix2, test)
matrix3 = as.matrix(matrix3)


#survival analysis
meta.survival$Gene = matrix3["JUNB",]
meta.survival$exp = "NA"
meta.survival[meta.survival$Gene >= mean(meta.survival$Gene),]$exp = "High"
meta.survival[meta.survival$Gene < mean(meta.survival$Gene),]$exp = "Low"
meta.survival$OS = as.numeric(meta.survival$OS)

fit = survfit(Surv(OS.time, OS)~exp, data = meta.survival)
#head(summary(fit))
a = ggsurvplot(fit, legend.title = "JUNB expression", palette = c("red","black"),
	pval=TRUE,
	risk.table = TRUE,
	tables.height = 0.2,
	tables.theme = theme_cleantable(),
	xlab = "Time(days)")

#gene correlation
cor.summary = matrix(NA, nrow=nrow(matrix3), ncol=2)
for (i in 1:nrow(matrix3)) {
	expr = matrix3[i,]
	tryCatch({
		a = cor.test(matrix3["SOX4",],expr)
		cor.summary[i,1] = a$estimate
		cor.summary[i,2] = a$p.value
		}, error = function(e){})
}
cor.summary = as.data.frame(cor.summary)
colnames(cor.summary) = c("beta","p")
rownames(cor.summary) = rownames(matrix3)
cor.summary = na.omit(cor.summary)
write.csv(cor.summary,file = "esca/cor_SOX4.csv")


#all gene survival
sur.summary = matrix(NA, nrow=nrow(matrix3), ncol=3)
for (i in 1:nrow(matrix3)) {
	meta.survival$exp = "NA"
	meta.survival$Gene = matrix3[i,]
	tryCatch({
		meta.survival[meta.survival$Gene >= mean(meta.survival$Gene),]$exp = "High"
		meta.survival[meta.survival$Gene < mean(meta.survival$Gene),]$exp = "Low"
		fit = survfit(Surv(OS.time, OS)~exp, data = meta.survival)
		fit.diff = survdiff(Surv(OS.time, OS)~exp, data = meta.survival)
		sur.summary[i,1:2] = fit$n
		sur.summary[i,3] = fit.diff$pvalue
		}, error = function(e){})
}

sur.summary = as.data.frame(sur.summary)
colnames(sur.summary) = c("N-expHigh","N-expLow","pvalue")
rownames(sur.summary) = rownames(matrix3)
sur.summary = na.omit(sur.summary)
write.csv(sur.summary, file = "esca/survival.scc.csv")

#cor plot 
library(ggpubr)
testdata = matrix3[c("ELK3","ETS1"),]
testdata2 = t(testdata)
testdata2 = as.data.frame(testdata2)
ggscatter(testdata2, x = "ETS1",y = "ELK3", add = "reg.line", palette = "jco") + 
	stat_cor(method = "pearson")


