##R
library(ggplot2)
library(ggpubr)

meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)
meta$myeloid_subtype = "NA"
meta[meta$myeloid_cluster == 0,]$myeloid_subtype = "Macro.C0"
meta[meta$myeloid_cluster == 1,]$myeloid_subtype = "Macro.C1"
meta[meta$myeloid_cluster == 2,]$myeloid_subtype = "cDC2.C2"
meta[meta$myeloid_cluster == 3,]$myeloid_subtype = "cDC2.C3"
meta[meta$myeloid_cluster == 4,]$myeloid_subtype = "Macro.C4"
meta[meta$myeloid_cluster == 5,]$myeloid_subtype = "Macro.C5"
meta[meta$myeloid_cluster == 6,]$myeloid_subtype = "Macro.C6"
meta[meta$myeloid_cluster == 7,]$myeloid_subtype = "Macro.C7"
meta[meta$myeloid_cluster == 8,]$myeloid_subtype = "Neutrophil.C8"
meta[meta$myeloid_cluster == 9,]$myeloid_subtype = "Macro.C9"
meta[meta$myeloid_cluster == 10,]$myeloid_subtype = "Macro.C10"
meta[meta$myeloid_cluster == 11,]$myeloid_subtype = "cDC1.C11"
meta[meta$myeloid_cluster == 12,]$myeloid_subtype = "Macro.C12"
meta[meta$myeloid_cluster == 13,]$myeloid_subtype = "Monocyte.C13"
meta[meta$myeloid_cluster == 14,]$myeloid_subtype = "Neutrophil.C14"
meta[meta$myeloid_cluster == 15,]$myeloid_subtype = "cDC3.C15"
meta[meta$myeloid_cluster == 16,]$myeloid_subtype = "Macro.C16"
meta[meta$myeloid_cluster == 17,]$myeloid_subtype = "cDC1.C17"
meta[meta$myeloid_cluster == 18,]$myeloid_subtype = "Macro.C18"
meta[meta$myeloid_cluster == 19,]$myeloid_subtype = "Macro.C19"
meta[meta$myeloid_cluster == 20,]$myeloid_subtype = "cDC2.C20"
meta[meta$myeloid_cluster == 21,]$myeloid_subtype = "Macro.C21"
meta[meta$myeloid_cluster == 22,]$myeloid_subtype = "Neutrophil.C22"
meta[meta$myeloid_cluster == 23,]$myeloid_subtype = "pDC.C23"
meta[meta$myeloid_cluster == 24,]$myeloid_subtype = "pDC.C24"
meta[meta$myeloid_cluster == 25,]$myeloid_subtype = "proliferating myeloid"
meta[meta$myeloid_cluster == 26,]$myeloid_subtype = "proliferating myeloid"

#1 tissue
test = table(x=meta$myeloid_subtype, y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("celltype","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

test1 = test[test$celltype %in% c("Monocyte.C13","Neutrophil.C8","Neutrophil.C14","Neutrophil.C22"),]
test2 = test[test$celltype %in% c("pDC.C24","cDC1.C11","cDC2.C3","cDC2.C20","cDC3.C15"),]
test3 = test[test$celltype %in% c("Macro.C21","Macro.C1","Macro.C5","Macro.C9","Macro.C6","Macro.C12","Macro.C16","Macro.C0"),]
test3$celltype = factor(test3$celltype, levels=c("Macro.C21","Macro.C1","Macro.C5","Macro.C0","Macro.C12","Macro.C9","Macro.C6","Macro.C16"))

p1 <- ggplot(data=test3,aes(x=celltype,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("figure/prop/macro.pdf",h=3,w=5)
p1
dev.off()

test2$celltype = factor(test2$celltype, levels=c("pDC.C24","cDC1.C11","cDC2.C3","cDC2.C20","cDC3.C15"))
p4 <- ggplot(data=test2,aes(x=tissue,y=Freq,fill=celltype))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#FF61C3","#00C083","#D39200","#BF80FF","#00B9E3"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("dc_cellabv.pdf",h=5,w=4)
p4
dev.off()

test3$celltype = factor(test3$celltype, levels=c("Macro.C21","Macro.C1","Macro.C5","Macro.C0","Macro.C12","Macro.C9","Macro.C6","Macro.C16"))
p3 <- ggplot(data=test3,aes(x=tissue,y=Freq,fill=celltype))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#DB72FB","#EE8045","#ACA300","#F8766D","#00C19F","#00BA38","#93AA00","#00B2F4"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("macro_cellabv.pdf",h=5,w=4)
p3
dev.off()


#2 patient
test = table(x=meta$myeloid_subtype, y=meta$patient)
test = as.data.frame(test)
colnames(test) = c("celltype","patient","Freq")
test$Freq = as.numeric(test$Freq)
test = test[test$celltype %in% c("Monocyte.C13","Neutrophil.C8","Neutrophil.C14","Neutrophil.C22","pDC.C24","cDC1.C11","cDC2.C3","cDC2.C20","cDC3.C15","Macro.C21","Macro.C1","Macro.C5","Macro.C9","Macro.C6","Macro.C12","Macro.C16","Macro.C0"),]
test$celltype = factor(test$celltype, levels=c("Monocyte.C13","Neutrophil.C8","Neutrophil.C14","Neutrophil.C22","pDC.C24","cDC1.C11","cDC2.C3","cDC2.C20","cDC3.C15","Macro.C21","Macro.C1","Macro.C5","Macro.C0","Macro.C12","Macro.C9","Macro.C6","Macro.C16"))

p2 <- ggplot(data=test,aes(x=celltype,y=Freq,fill=patient))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#008000","#BC80BD"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#3
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$myeloid_cluster == 16,]
test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref[rownames(test),]
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)
test1$tissue = factor(test1$tissue, levels = c("Tumor","Adjacent","Distant"))

my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/prop/C16.pdf",h=5,w=4)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()



#4 smoker
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$myeloid_celltype == "cDC1",]
meta1 = meta[meta$myeloid_celltype == "cDC2",]
meta1 = meta[meta$myeloid_celltype == "cDC3",]
meta1 = meta[meta$myeloid_celltype == "M1-like Mac",]
meta1 = meta[meta$myeloid_celltype == "M2-like Mac",]
meta1 = meta[meta$myeloid_celltype == "Monocyte",]
meta1 = meta[meta$myeloid_celltype == "Neutrophil",]



test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)


test1$tissue = c("adj-smo","adj-smo","adj-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-smo","adj-non-smo","adj-smo","adj-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-smo","dist-smo","dist-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-smo","tumor-smo","tumor-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo")
test1$tissue = factor(test1$tissue, levels = c("tumor-smo","tumor-non-smo","adj-smo","adj-non-smo","dist-smo","dist-non-smo"))
my_comparison = list(c("tumor-smo","tumor-non-smo"),c("adj-smo","adj-non-smo"),c("dist-smo","dist-non-smo"))


pdf("figure/prop-smk-cdc3.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()


test1$tissue = c("adj-drk","adj-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk")
test1$tissue = factor(test1$tissue, levels = c("tumor-drk","tumor-non-drk","adj-drk","adj-non-drk","dist-drk","dist-non-drk"))
my_comparison = list(c("tumor-drk","tumor-non-drk"),c("adj-drk","adj-non-drk"),c("dist-drk","dist-non-drk"))


pdf("figure/prop-drk-cdc3.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()
