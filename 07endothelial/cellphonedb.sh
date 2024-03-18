###准备流程
##1 安装anaconda
wget -c https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
bash Anaconda2-2019.03-Linux-x86_64.sh
#添加到环境变量
export PATH="/share/home/daijiacheng/anaconda2/bin:$PATH"
#添加bioconda channels
conda config --add channels bioconda

##2 配置cellphoneDB环境
conda create -n cellphonedb3 python=3.7
source activate cellphonedb3
pip install cellphonedb

##3 提取seurat对象
library(Seurat)

data1 = readRDS("../02EC_cluster/new0423/02ECdataset.16580.rds")
data2 = readRDS("../08fibroblast/08FB.sct.rds")
data3 = readRDS("../08smc/08SMC.sct.rds") 

genelist = Reduce(intersect, list(rownames(data1),rownames(data2),rownames(data3)))
#17980

meta1 = data1@meta.data[data1@meta.data$celltype %in% c("Artery_UNC5B","tip cell_COL4A1+","stalk cell_RGCC+","Vein_SELE+","Vein_ISG15+"),c("celltype","seurat_clusters")]
meta2 = data2@meta.data[data2@meta.data$celltype %in% c("CAF_1","CAF_2","CAF_3","CAF_4","CAF_5","CAF_6"),c("celltype","seurat_clusters")]
meta3 = data3@meta.data[data3@meta.data$celltype %in% c("Pericyte_1","Pericyte_2","SMC_8"),c("celltype","seurat_clusters")]
matrix1 = data1@assays$SCT@data[genelist,rownames(meta1)]
matrix2 = data2@assays$SCT@data[genelist,rownames(meta2)]
matrix3 = data3@assays$SCT@data[genelist,rownames(meta3)]
data.input = cbind(matrix1, matrix2, matrix3)
gene = rownames(data.input)
rt1 = cbind(gene, as.matrix(data.input))
write.table(rt1, "07cellphonedb/cellphonedb_count.txt", sep = "\t", quote = F, row.names = F, col.names = T)
meta = rbind(meta1, meta2, meta3)
meta[,2] = meta[,1]
meta[,1] = rownames(meta)
colnames(meta) = c("sampleID","celltype")
write.table(meta, "07cellphonedb/cellphonedb_meta.txt", sep = "\t", quote = F, row.names = F)


##4 cellphoneDB
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads 32 --counts-data gene_name
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot cellphonedb_meta.txt

#plot
#VEGF
library(tidyverse)
mypvals = read.delim("pvalues.txt", check.names=FALSE)
mymeans = read.delim("means.txt", check.names=FALSE)

VEGF = grep("VEGF",mymeans$interacting_pair, value=T)


mymeans %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair",starts_with("Fibroblast"),ends_with("Fibroblast")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair",starts_with("Fibroblast"),ends_with("Fibroblast")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")

pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")

summary((filter(pldf,means>0.1))$means)

pldf %>% filter(means > 0.1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


#ACKR1
ACKR1 = grep("ACKR1",mymeans$interacting_pair, value=T)

mymeans %>% dplyr::filter(interacting_pair %in% ACKR1) %>% dplyr::select("interacting_pair",starts_with("T cell"),ends_with("T cell")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% ACKR1) %>% dplyr::select("interacting_pair",starts_with("T cell"),ends_with("T cell")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")

pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")

summary((filter(pldf,means>0.1))$means)

pldf %>% filter(means > 0.1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))