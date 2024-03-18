##R 
##single cell data of ESCC
##step4.5: SCENIC
##Jiacheng Dai #daicy0424@gmail.com
##2022/09/09

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
rm(list=ls())

setwd("01ECatlas/06SCENIC")
dir.create("int")
data = readRDS("../02EC_cluster/02ECdataset.16580.rds")
meta=read.csv("total.16580.csv",header=T,row.names=1)
cellinfo = meta[,c("tissue","celltype")]
colnames(cellinfo) = c("tissue","celltype")
saveRDS(cellinfo, file="int/cellInfo.rds")

#Initialize settings
cellinfo = cellinfo[cellinfo$tissue=="Tumor",]
exprMat = as.matrix(data@assays$RNA@counts[,rownames(cellinfo)])
mydbDIR = "./hg38"
mydbs = c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) = c("10kb","500bp")
scenicOptions = initializeScenic(org="hgnc",
                                nCores=16,
                                dbDir=mydbDIR,
                                dbs=mydbs,
                                datasetTitle="combine")
saveRDS(scenicOptions, "int/scenicOp-EC.rds")

#co-expression network
genesKept = geneFiltering(exprMat, scenicOptions, minCountsPerGene = 3*0.01*ncol(exprMat), minSamples=ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

saveRDS(cellinfo, file="int/cellInfo.rds")
saveRDS(scenicOptions, "int/scenicOp-EC.rds")

#build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

#binarize activity
#aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
#savedSelections <- shiny::runApp(aucellApp)
#newThresholds <- savedSelections$thresholds
#scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

#export
saveRDS(cellinfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)
saveRDS(cellinfo, file="int/cellInfo.rds")
saveRDS(scenicOptions, "int/scenicOp-EC.rds")

##visualize
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(aucell_regulonAUC),cellAnnotation=cellinfo[colnames(aucell_regulonAUC),"celltype"])
rssPlot <- plotRSS(rss)
rssPlot$plot
pdf("rssplot.pdf")
rssPlot$plot
dev.off()

#Complex Heatmap
regulonAUC=loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC=regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType = sapply(split(rownames(cellinfo),cellinfo$celltype),function(cells)rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pdf("int/ComplexHeatmap.pdf")
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()


# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("ATF3", "JUND", "FOSB","JUNB", "JUN","FOS","STAT1","IRF7")],], plots="Expression")
pdf("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()


#画scenic的那些regulon的表达情况
pdf("scenic.vlnplot/irf7.pdf",h=5,w=8)
VlnPlot(data, features="IRF7", pt.size=0)
dev.off()
