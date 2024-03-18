#R
library(ggplot2)

meta = read.csv("metadata.res1.3.csv",header=T,row.names=1)
#1 tissue
test = table(x=meta$epi_cluster,y=meta$tissue)
test = as.data.frame(test)
colnames(test) = c("cluster","tissue","Freq")
test$Freq = as.numeric(test$Freq)
test$tissue = factor(test$tissue, levels = c("Tumor","Adjacent","Distant"))

tumor = test[test$cluster %in% c(1,2,3,4,5,6,7,8,11,12,13,14,15,16,18,20,21,25),]
dist = test[test$cluster %in% c(0,9,10,17,19,22,23,24),]

p1 <- ggplot(data=dist,aes(x=cluster,y=Freq,fill=tissue))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#ff0000","#008000","#0000FF"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

#2 patient
test = table(x=meta$epi_cluster,y=meta$patient)
test = as.data.frame(test)
colnames(test) = c("cluster","patient","Freq")
test$Freq = as.numeric(test$Freq)

tumor = test[test$cluster %in% c(1,2,3,4,5,6,7,8,11,12,13,14,15,16,18,20,21,25),]
dist = test[test$cluster %in% c(0,9,10,17,19,22,23,24),]

p2 <- ggplot(data=dist,aes(x=cluster,y=Freq,fill=patient))+
     geom_bar(stat="identity",position="fill")+
     theme_minimal()+
     scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#008000","#BC80BD"))+
     theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


#3 
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






