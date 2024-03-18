#R
library(ggplot2)
library(ggpubr)

meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)
#1 tissue
test = table(x=meta$fb_celltype, y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("celltype","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

p1 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("figure/cellprop.pdf",h=4,w=11)
p1
dev.off()

p2 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("figure/cellvalue.pdf",h=4,w=11)
p2
dev.off()

#2
test = table(x=meta$fb_celltype, y=meta$patient)
test = as.data.frame(test)
colnames(test) = c("celltype","patient","Freq")
test$Freq = as.numeric(test$Freq)

p2 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=patient))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#008000","#BC80BD"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("figure/patient_prop.pdf",h=6,w=10)
p2
dev.off()


#3 boxplot
ref = table(x=meta$Patient_ID, y=meta$tissue)

meta1 = meta[meta$Fib_celltype == "Fib18",]
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


#4
meta[meta$fb_cluster %in% c(6,14,16,18,19,25,26,27),]$major = "CAF"
test = table(x=meta$major,y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("cluster","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

p1 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
