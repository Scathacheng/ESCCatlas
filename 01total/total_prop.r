#R
library(ggplot2)
library(ggpubr)

meta = read.csv("meta_allsample.csv",header=T,row.names=1)
#1 tissue
test = table(x=meta$major, y=meta$tissue)
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

#2 patient
test = table(x=meta$patient, y=meta$tissue)
test = as.data.frame(test)
colnames(test)=c("patient","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

p1 <- ggplot(data=test,aes(x=patient,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="dodge")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p2 <- ggplot(data=test,aes(x=patient,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#3
ref = table(x=meta$patient, y=meta$tissue)
meta1 = meta[meta$major=="T cell",]
test = table(x=meta1$patient, y=meta1$tissue)
test1 = test/ref
test1 = as.data.frame(test1)
test1 = na.omit(test1)
colnames(test1) = c("patient","tissue","Proportion")
test1$Proportion = as.numeric(test1$Proportion)
test1$tissue = factor(test1$tissue, levels = c("Tumor","Adjacent","Distant"))

my_comparison = list(c("Tumor","Adjacent"),c("Tumor","Distant"),c("Adjacent","Distant"))

ggboxplot(test1, x="tissue",y="Proportion", color="tissue",add="jitter")+
	stat_compare_means(comparisons = my_comparison)


