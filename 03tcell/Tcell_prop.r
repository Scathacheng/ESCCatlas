#R
library(ggplot2)
library(ggpubr)

meta = read.csv("metadata.res1.4.csv",header=T,row.names=1)
#1 tissue
test = table(x=meta$Tcellsubtype, y=meta$tissue)
rownames(test) = c("naive CD4","naive CD8","exhausted T cell","naive T cell","cytotoxic T cell","proliferating T cell","DNAJB1 T cell","IGFBP7 T cell","CXCR4 Tem","GZMK Tem","Treg")
test = as.data.frame(test)
colnames(test) = c("cluster","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

cd4tcell = test[test$cluster %in% c(1,2,14,18),]
cd8tcell = test[test$cluster %in% c(0,4,6,9,12,13,15),]
other = test[test$cluster %in% c(3,5,7,8,10,11,16,17,19),]

p1 <- ggplot(data=test,aes(x=cluster,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=cluster,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

pdf("cellabv.pdf",h=3,w=6)
p2
dev.off()

#"naive CD4","naive CD8","exhaust T","naive T","cytotoxic T","proliferating T","DNAJB1 T","IGFBP7 T","CXCR4 Tem","GZMK Tem","Treg"
#"#00BADE","#00A6FF","#64B200","#B385FF","#DB8E00","#EF67EB","#AEA200","#00C1A7","#F8766D","#00BD5C","#FF63B6"

test$cluster = factor(test$cluster, levels=c("naive T cell","naive CD4","Treg","proliferating T cell","DNAJB1 T cell","IGFBP7 T cell","naive CD8","GZMK Tem","CXCR4 Tem","cytotoxic T cell","exhausted T cell"))

p3 <- ggplot(data=test,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#B385FF","#00BADE","#FF63B6","#EF67EB","#AEA200","#00C1A7","#00A6FF","#00BD5C","#F8766D","#DB8E00","#64B200"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#h=6,w=3
pdf("tissue_cellabv.pdf",h=6,w=5)
p3
dev.off()

exhaust = test[test$cluster %in% c("exhausted T cell","cytotoxic T cell","CXCR4 Tem","GZMK Tem"),]
exhaust$cluster = factor(exhaust$cluster, levels=c("GZMK Tem","CXCR4 Tem","cytotoxic T cell","exhausted T cell"))
p4 <- ggplot(data=exhaust,aes(x=tissue,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#00BD5C","#F8766D","#DB8E00","#64B200"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pdf("exhaust_cellprop.pdf",h=6,w=5)
p4
dev.off()

#3
test = table(x=meta$Tcellsubtype, y=meta$patient)
rownames(test) = c("naive CD4","naive CD8","exhausted T cell","naive T cell","cytotoxic T cell","proliferating T cell","DNAJB1 T cell","IGFBP7 T cell","CXCR4 Tem","GZMK Tem","Treg")
test = as.data.frame(test)
colnames(test) = c("cluster","patient","Number")
test$Number = as.numeric(test$Number)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

p2 <- ggplot(data=test,aes(x=patient,y=Number,fill=cluster))+
     geom_bar(stat="identity",position="stack")+
     theme_minimal()+
     scale_fill_manual(values=c("#00BADE","#00A6FF","#64B200","#B385FF","#DB8E00","#EF67EB","#AEA200","#00C1A7","#F8766D","#00BD5C","#FF63B6"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#h=6,w=8

#2 
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$tcell_cluster==6,]
test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref[rownames(test),]
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)
test1$tissue = factor(test1$tissue, levels = c("Tumor","Adjacent","Distant"))
my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

pdf("figure/C6.pdf",h=5,w=4)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()



#
p1 <- ggplot(data=cd8tcell,aes(x=patient,y=Freq,fill=cluster))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#1F8B43","#252F6B","#D72029","#8AC860","#6D68AD","#F38AA6","#FCB83D"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


#smoker
ref = table(x=meta$patient, y=meta$tissue)

meta1 = meta[meta$tcell_celltype %in% c("Treg.CCR6.C14","Treg.CXCR3.C22","Treg.TNFRSF4.C20"),]
meta1 = meta[meta$tcell_celltype %in% c("Tem.CXCR4.C19","Tem.GZMK.C16","Tem.GZMK.C21","Tem.GZMK.C9","Tem.GZMK.C28"),]
meta1 = meta[meta$tcell_celltype %in% c("TNK.C24","TNK.C27","TNK.C4"),]
meta1 = meta[meta$tcell_celltype %in% c("Tex.C12","Tex.C8"),]
meta1 = meta[meta$tcell_celltype %in% c("Tn.CD8.C0","Tn.CD8.C1","Tn.CD8.C18","Tn.CD8.C25"),]
meta1 = meta[meta$tcell_celltype %in% c("Tcell.DNAJB1.C10","Tcell.DNAJB1.C29","Tcell.DNAJB1.C5"),]
meta1 = meta[meta$tcell_celltype %in% c("Tcell.IGFBP7.C13","Tcell.IGFBP7.C23"),]
meta1 = meta[meta$tcell_celltype == "Tn.CD4.C6",]



test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)


test1$tissue = c("adj-smo","adj-smo","adj-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-non-smo","adj-smo","adj-non-smo","adj-smo","adj-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-smo","dist-smo","dist-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-non-smo","dist-smo","dist-non-smo","dist-smo","dist-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-smo","tumor-smo","tumor-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-non-smo","tumor-smo","tumor-non-smo","tumor-smo","tumor-non-smo")
test1$tissue = factor(test1$tissue, levels = c("tumor-smo","tumor-non-smo","adj-smo","adj-non-smo","dist-smo","dist-non-smo"))

my_comparison = list(c("tumor-smo","tumor-non-smo"),c("adj-smo","adj-non-smo"),c("dist-smo","dist-non-smo"))
pdf("figure/prop-smk-tnCD4.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()

test1$tissue = c("adj-drk","adj-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-non-drk","adj-drk","adj-non-drk","adj-non-drk","adj-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-non-drk","dist-drk","dist-non-drk","dist-non-drk","dist-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk","tumor-drk","tumor-non-drk","tumor-non-drk","tumor-non-drk")
test1$tissue = factor(test1$tissue, levels = c("tumor-drk","tumor-non-drk","adj-drk","adj-non-drk","dist-drk","dist-non-drk"))

my_comparison = list(c("tumor-drk","tumor-non-drk"),c("adj-drk","adj-non-drk"),c("dist-drk","dist-non-drk"))
pdf("figure/prop-drk-tnCD4.pdf",h=5,w=8)
ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)
dev.off()
