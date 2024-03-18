#R
#EC 
library(ggplot2)
library(ggpubr)

meta=read.csv("total.16580.csv",header=T,row.names=1)

#boxplot
ref = table(x=meta$Patient_ID, y=meta$tissue)
table(meta$celltype)
meta1 = meta[meta$celltype == "Vein_SELE+",]
test = table(x=meta1$Patient_ID, y=meta1$tissue)
test1 = test/ref[rownames(test),]
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)
test1$tissue = factor(test1$tissue, levels = c("Tumor","Adjacent","Distant"))

my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)

#
test = table(x=meta$celltype, y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("celltype","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

p1 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


#4
#NSCLC
meta=read.csv("Metadata.csv",header=T,row.names=1)
test = table(x=meta$Cluster,y=meta$NEC.TEC)
test = as.data.frame(test)
colnames(test) = c("cluster","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("TEC","NEC"))

p1 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#0000FF","#B385FF","#00BADE","#FF63B6","#EF67EB","#AEA200","#00C1A7","#00A6FF","#00BD5C","#F8766D","#DB8E00","#64B200","#ff0000"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#0000FF","#B385FF","#00BADE","#FF63B6","#EF67EB","#AEA200","#00C1A7","#00A6FF","#00BD5C","#F8766D","#DB8E00","#64B200","#ff0000"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#ESCC
meta=read.csv("total.16580.csv",header=T,row.names=1)
test = table(x=meta$celltype,y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("cluster","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))
test$cluster = factor(test$cluster, levels = c("Artery_UNC5B","Artery_DEPP1","Artery_GAS6","Capillary_CA4","Capillary_FABP5","Capillary_SOCS3","EC_HSP","EC_IGKC","PCV_ACTG1","Vein_CPE","Vein_HLA","Vein_IL6","Vein_ISG15","Vein_SELE","stalk cell_RGCC","tip cell_COL4A1"))

p1 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))