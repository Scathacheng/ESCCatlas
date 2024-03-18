library(GSVA)
library(pheatmap)


#注释32个上皮细胞亚群
###1 load the cell marker list
#1 nature communication
genelist = read.csv("cellident2.csv",header=T)
colnames(genelist) = c("Celltype","Gene")
list = split(as.matrix(genelist)[,2], genelist[,1])

#2 cancer cell
genelist = read.csv("cellident1.csv",header=T)
colnames(genelist) = c("Celltype","Gene")
list = split(as.matrix(genelist)[,2], genelist[,1])

###2 load our data
#已经计算过cluster均值的基因表达矩阵
data = read.csv("expr.csv",header=T,row.names=1)
colnames(data) = paste0('C',0:31)

tumor = data[,c("C0","C1","C2","C3","C6","C8","C9","C10","C11","C12","C14","C16","C17","C18","C19","C20","C23","C25","C27","C28","C29","C30","C31")]
dist = data[,c("C4","C5","C7","C13","C15","C21","C22","C24","C26")]

data.gsva = gsva(as.matrix(dist), list, method="ssgsea", kcdf = "Gaussian", abs.ranking=TRUE)
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

