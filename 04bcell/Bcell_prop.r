#R
library(ggplot2)
library(ggpubr)

meta = read.csv("metadata.res1.1.csv",header=T,row.names=1)
#1 tissue
test = table(x=meta$Bcell_identity, y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("celltype","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

iga = test[test$celltype %in% c("IgA PC1","IgA PC2","IgA PC3","IgA PC4","IgA PC5","IgA PC6","IgA PC7","IgA PC8","IgM PC"),]
igg = test[test$celltype %in% c("IgG3 PC1","IgG3 PC2","IgG3 PC3","IgG3 PC4","IgG3 PC5","IgG3 PC6","IgG3 PC7","IgG3 PC8","IgG3 PC9","IgG2 PC1","IgG2 PC2","IgG2 PC3","IgG2 PC4","IgG4 PC1","IgG4 PC2"),]
bcell = test[test$celltype %in% c("Naive B cell","Mem B1","Mem B2","Mem B3","Mem B4","Mem B5","Mem B6","Mem B7","Mem B8","Stress B1","Stress B2"),]

p1 <- ggplot(data=bcell,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
p2 <- ggplot(data=bcell,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))



pdf("figure/cellprop_bcell.pdf",h=3,w=6)
p1
dev.off()
pdf("figure/cellabv_bcell.pdf",h=3,w=6)
p2
dev.off()


#3 
test = table(x=meta$bcell_celltype, y=meta$patient)
test = as.data.frame(test)
colnames(test) = c("celltype","patient","Freq")
test$Freq = as.numeric(test$Freq)

iga = test[test$celltype %in% c("IgA PC1","IgA PC2","IgA PC3","IgA PC4","IgA PC5","IgA PC6","IgA PC7","IgA PC8","IgM PC"),]
igg = test[test$celltype %in% c("IgG3 PC1","IgG3 PC2","IgG3 PC3","IgG3 PC4","IgG3 PC5","IgG3 PC6","IgG3 PC7","IgG3 PC8","IgG3 PC9","IgG2 PC1","IgG2 PC2","IgG2 PC3","IgG2 PC4","IgG4 PC1","IgG4 PC2"),]
bcell = test[test$celltype %in% c("Naive B cell","Mem B1","Mem B2","Mem B3","Mem B4","Mem B5","Mem B6","Mem B7","Mem B8","Stress B1","Stress B2"),]

p2 <- ggplot(data=bcell,aes(x=celltype,y=Freq,fill=patient))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#008000","#BC80BD"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("figure/patientprop_bcell.pdf",h=6,w=6)
p2
dev.off()


#2 boxplot
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$bcell_celltype == "IgG4 PC2",]
test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref[rownames(test),]
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)
test1$tissue = factor(test1$tissue, levels = c("Tumor","Adjacent","Distant"))

my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/ig4pc2.pdf",h=5,w=4)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()



#4 smoker
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$bcell_celltype %in% c("IgG3 PC1","IgG3 PC2","IgG3 PC3","IgG3 PC4","IgG3 PC5","IgG3 PC6","IgG3 PC7","IgG3 PC8","IgG3 PC9","IgG2 PC1","IgG2 PC2","IgG2 PC3","IgG2 PC4","IgG4 PC1","IgG4 PC2"),]
meta1 = meta[meta$bcell_celltype %in% c("IgA PC1","IgA PC2","IgA PC3","IgA PC4","IgA PC5","IgA PC6","IgA PC7","IgA PC8","IgM PC"),]
meta1 = meta[meta$bcell_celltype %in% c("Naive B cell","Mem B1","Mem B2","Mem B3","Mem B4","Mem B5","Mem B6","Mem B7","Mem B8","Stress B1","Stress B2"),]


test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)


test1$tissue = c("adj-smo","adj-smo","adj-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-smo","adj-non-smo","adj-smo","adj-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-smo","dist-smo","dist-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-smo","tumor-smo","tumor-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo")
test1$tissue = factor(test1$tissue, levels = c("tumor-smo","tumor-non-smo","adj-smo","adj-non-smo","dist-smo","dist-non-smo"))

my_comparison = list(c("tumor-smo","tumor-non-smo"),c("adj-smo","adj-non-smo"),c("dist-smo","dist-non-smo"))
pdf("figure/prop-smk-bcell.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()

test1$tissue = c("adj-drk","adj-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk")
test1$tissue = factor(test1$tissue, levels = c("tumor-drk","tumor-non-drk","adj-drk","adj-non-drk","dist-drk","dist-non-drk"))

my_comparison = list(c("tumor-drk","tumor-non-drk"),c("adj-drk","adj-non-drk"),c("dist-drk","dist-non-drk"))
pdf("figure/prop-drk-bcell.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()
