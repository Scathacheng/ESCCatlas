#R
library(tidyverse)
library(dplyr)
library(reshape2)


####chord plot
#读取
tumor = read.delim("tumor/significant_means.txt",check.names=FALSE)
adj = read.delim("adj/significant_means.txt",check.names=FALSE)
dist = read.delim("dist/significant_means.txt",check.names=FALSE)

test1 = tumor[tumor$rank < 0.01,] 
test2 = adj[adj$rank < 0.01,] 
test3 = dist[dist$rank < 0.01,] 


test4 = test1[,colSums(is.na(test1)) != nrow(test1)] # 314 9985
test5 = test2[,colSums(is.na(test2)) != nrow(test2)] # 319 8568
test6 = test3[,colSums(is.na(test3)) != nrow(test3)] # 344 8200

mytumor = test4[,c(-1,-3:-12)]
myadj = test5[,c(-1,-3:-12)]
mydist = test6[,c(-1,-3:-12)]

mytumor = na.omit(melt(mytumor))
myadj = na.omit(melt(myadj))
mydist = na.omit(melt(mydist))

colnames(mytumor) = colnames(myadj) = colnames(mydist) = c("interacting_pair","Cell_Cell","means")

write.csv(mytumor, "chord/mytumor.csv")
write.csv(myadj, "chord/myadj.csv")
write.csv(mydist, "chord/mydist.csv")

library(circlize)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(reshape2)


#去重和分列
mytumor = read.csv("chord/mytumor.csv",header=T,row.names=1)
result1 = table(mytumor$Cell_Cell)
result1 = as.data.frame(result1)
write.csv(result1, "chord/chord_tumor.csv")

myadj = read.csv("chord/myadj.csv",header=T,row.names=1)
result2 = table(myadj$Cell_Cell)
result2 = as.data.frame(result2)
write.csv(result2,"chord/chord_adj.csv")

mydist = read.csv("chord/mydist.csv",header=T,row.names=1)
result3 = table(mydist$Cell_Cell)
result3 = as.data.frame(result3)
write.csv(result3,"chord/chord_dist.csv")

#chord plot
tumor = read.csv("chord/chord_tumor.csv",header=T,row.names=1)

test = dcast(tumor, start~end)
rownames(test) = test[,1]
test = test[,-1]
set.seed(888)
grid.col <- brewer.pal(4, "Set1")
chordDiagram(test, grid.col = grid.col,
    annotationTrack = c("name","grid"),
    annotationTrackHeight = c(0.03, 0.01))



adj = read.csv("chord/chord_adj.csv",header=T,row.names=1)
test = dcast(adj, start~end)
rownames(test) = test[,1]
test = test[,-1]
set.seed(888)
grid.col <- brewer.pal(4, "Set1")
chordDiagram(test, grid.col = grid.col,
    annotationTrack = c("name","grid"),
    annotationTrackHeight = c(0.03, 0.01))


dist = read.csv("chord/chord_dist.csv",header=T,row.names=1)
test = dcast(dist, start~end)
rownames(test) = test[,1]
test = test[,-1]
set.seed(888)
grid.col <- brewer.pal(4, "Set1")
chordDiagram(test, grid.col = grid.col,
    annotationTrack = c("name","grid"),
    annotationTrackHeight = c(1, 1))



####heatmap
tumor = read.delim("tumor/significant_means.txt",check.names=FALSE)
adj = read.delim("adjacent/significant_means.txt",check.names=FALSE)
dist = read.delim("distant/significant_means.txt",check.names=FALSE)

test1 = tumor[tumor$rank < 0.01,] 
test2 = adj[adj$rank < 0.01,] 
test3 = dist[dist$rank < 0.01,] 

test4 = test1[,colSums(is.na(test1)) != nrow(test1)] #342 7337
test5 = test2[,colSums(is.na(test2)) != nrow(test2)] #308 6185
test6 = test3[,colSums(is.na(test3)) != nrow(test3)] #319 6029

mytumor = test4[,c(-1,-3:-12)]
myadj = test5[,c(-1,-3:-12)]
mydist = test6[,c(-1,-3:-12)]

mytumor = na.omit(melt(mytumor))
myadj = na.omit(melt(myadj))
mydist = na.omit(melt(mydist))

colnames(mytumor) = colnames(myadj) = colnames(mydist) = c("interacting_pair","Cell_Cell","means")


#tumor specific interaction
mytumor = read.csv("mytumor.csv",header=T,row.names=1)
myadj = read.csv("myadj.csv",header=T,row.names=1)
mydist = read.csv("mydist.csv",header=T,row.names=1)
specific = setdiff(setdiff(mytumor$interacting_pair, myadj$interacting_pair),mydist$interacting_pair)
test = mytumor[mytumor$interacting_pair %in% specific,]

#1 vascular 
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = c("KDR_VEGFC","PDGFRB_PDGFD")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", ends_with("+")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", ends_with("+")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.5) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#2 FGF
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = c("FGF2_FGFR2","FGFR1_FGF3")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Fib","CAF"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Fib","CAF"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.56) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#3 cDC
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = c("ACKR1_CCL17","CCR4_CCL17")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Treg","Vein")), ends_with(c("cDC2.C3","cDC3.C15"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Treg","Vein")), ends_with(c("cDC2.C3","cDC3.C15"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 1.5) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#4 Treg
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = c("CCR4_CCL17","CCL22_CCR4")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Treg.TNFRSF4.C20","Treg.CXCR3.C22","Treg.CCR6.C14","cDC3.C15","cDC2.C3")), ends_with(c("Treg.TNFRSF4.C20","Treg.CXCR3.C22","Treg.CCR6.C14","cDC3.C15","cDC2.C3"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Treg.TNFRSF4.C20","Treg.CXCR3.C22","Treg.CCR6.C14","cDC3.C15","cDC2.C3")), ends_with(c("Treg.TNFRSF4.C20","Treg.CXCR3.C22","Treg.CCR6.C14","cDC3.C15","cDC2.C3"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.1810) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#4 SPP1
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)

mymeans %>% dplyr::filter(interacting_pair == "SPP1_CCR8") %>% dplyr::select("interacting_pair", starts_with("Macro"), ends_with(c("cDC3.C15","cDC2.C3"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair == "SPP1_CCR8") %>% dplyr::select("interacting_pair", starts_with("Macro"), ends_with(c("cDC3.C15","cDC2.C3"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 1.3) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


#normal specific 
#1 
rm(list=ls())
mypvals = read.delim("adj/pvalues.txt", check.names=FALSE)
mymeans = read.delim("adj/means.txt", check.names=FALSE)
genelist = c("PDGFA_PDGFRA","NOTCH3_JAG2")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair",starts_with("F"),ends_with("+")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair",starts_with("F"),ends_with("+")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
#adj >0.5
#dist >0.57
pldf %>% filter(means > 0.57) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#2
rm(list=ls())
mypvals = read.delim("dist/pvalues.txt", check.names=FALSE)
mymeans = read.delim("dist/means.txt", check.names=FALSE)
genelist = c("LTA_TNFRSF1A","FLT3_FLT3LG")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair",starts_with(c("cDC1","cDC2","Macro"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair",starts_with(c("cDC1","cDC2","Macro"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
#adj > 0.25
#dist > 0.25
pldf %>% filter(means > 0.25) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


#3
rm(list=ls())
mypvals = read.delim("adj/pvalues.txt", check.names=FALSE)
mymeans = read.delim("adj/means.txt", check.names=FALSE)
genelist = c("FCGR2A_CXCL9")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("cDC1","cDC2","cDC3"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("cDC1","cDC2","cDC3"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
#adj > 0.5
#dist > 0.47
pldf %>% filter(means > 0.5) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))



#angiogenesis
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = grep("TGF",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", ends_with(c("HLA+","CPE+","IL6+","SELE+"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", ends_with(c("HLA+","CPE+","IL6+","SELE+"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.3) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


#VEGF tumor-adj-dist=0.7
#PDGF 0.5
#TGF 0.3

#exhausted T cell
rm(list=ls())
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = grep("CYSLTR1",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#CXCL13 1.9
#TNFRSF4 0.8
#PDCD1 0.4
#TIGIT 0.4
#CTLA


genelist = grep("CTLA",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.3) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#T cell
rm(list=ls())
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = grep("CYSLTR1",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with(c("Tex.C12","Tex.C8"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


#MMP
#TGF 很多
#BMP
#IL
#TNF
#INF
#CXCL
#RTK

rm(list=ls())
mypvals = read.delim("dist/pvalues.txt", check.names=FALSE)
mymeans = read.delim("dist/means.txt", check.names=FALSE)
genelist = grep("TNF",mymeans$interacting_pair, value=T)
genelist = c("TNFRSF1B_GRN","TNFRSF1A_GRN")
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with("Epi")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with("Epi")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


rm(list=ls())
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genelist = grep("FGF",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with("CAF")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genelist) %>% dplyr::select("interacting_pair", starts_with("CAF")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary(pldf$means)
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 1) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))







#example
table(test$interacting_pair)
#VEGF
rm(list=ls())
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
VEGF = c("VEGFA_GRIN2B","KDR_VEGFC","EPHB6_EFNB2","PDGFR complex_PDGFD")
mymeans %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair", ends_with("tip cell_COL4A1+")) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair", ends_with("tip cell_COL4A1+")) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary((filter(pldf,means>0.1))$means)
pldf %>% filter(means > 0.21) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

#DC
rm(list=ls())
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
ackr1 = c("ACKR1_CCL17")
mymeans %>% dplyr::filter(interacting_pair %in% ackr1) %>% dplyr::select("interacting_pair", ends_with(c("cDC2.C3","cDC3.C15"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% ackr1) %>% dplyr::select("interacting_pair", ends_with(c("cDC2.C3","cDC3.C15"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary((filter(pldf,means>0.2))$means)
pldf %>% filter(means > 0.2) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))
#Tex
mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
genes = names(table(test$interacting_pair))
mymeans %>% dplyr::filter(interacting_pair %in% genes) %>% dplyr::select("interacting_pair", starts_with("Tex"), ends_with(c("Tex.C8","Tex.C12"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% genes) %>% dplyr::select("interacting_pair", starts_with("Tex"),ends_with(c("Tex.C8","Tex.C12"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")
pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")
summary((filter(pldf,means>0.2))$means)
pldf %>% filter(means > 0.5) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))


mypvals = read.delim("tumor/pvalues.txt", check.names=FALSE)
mymeans = read.delim("tumor/means.txt", check.names=FALSE)
VEGF = grep("VEGF",mymeans$interacting_pair, value=T)
mymeans %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair", ends_with(c("COL4A1+","RGCC+","UNC5B+"))) %>% reshape2::melt() -> meansdf
colnames(meansdf) <- c("interacting_pair","Cell_Cell","means")
mypvals %>% dplyr::filter(interacting_pair %in% VEGF) %>% dplyr::select("interacting_pair", ends_with(c("COL4A1+","RGCC+","UNC5B+"))) %>% reshape2::melt() -> pvalsdf
colnames(pvalsdf) <- c("interacting_pair","Cell_Cell","pvals")

pvalsdf$joinlab = paste0(pvalsdf$interacting_pair,"_",pvalsdf$Cell_Cell)
meansdf$joinlab = paste0(meansdf$interacting_pair,"_",meansdf$Cell_Cell)
pldf = merge(pvalsdf, meansdf, by="joinlab")

summary((filter(pldf,means>0.1))$means)

pldf %>% filter(means > 0.7) %>% 
        ggplot(aes(Cell_Cell.x, interacting_pair.x)) + 
        geom_point(aes(color=means, size=-log10(pvals+0.0001))) + 
        scale_size_continuous(range=c(1,3)) + 
        scale_color_gradient2(high="red", mid="yellow", low="darkblue", midpoint=25) +
        theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust=-0.1,vjust=0.8))

