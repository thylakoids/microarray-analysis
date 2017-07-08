rm(list=ls())
library(topGO)
library("org.Hs.eg.db")
library(ggplot2)
goAnn <- get("org.Hs.egGO")
universe <- Lkeys(goAnn)
#genes <- c("AREG", "FKBP5", "CXCL13", "KLF9", "ZC3H12A", "P4HA1", "TLE1", "CREB3L2", "TXNIP", "PBX1", "GJA1", "ITGB8", "CCL3", "CCND2", "KCNJ15", "CFLAR", "CXCL10", "CYSLTR1", "IGFBP7", "RHOB", "MAP3K5", "CAV2", "CAPN2", "AKAP13", "RND3", "IL6ST", "RGS1", "IRF4", "G3BP1", "SEL1L", "VEGFA", "SMAD1", "CCND1", "CLEC3B", "NEB", "AMD1", "PDCD4", "SCD", "TM2D3", "BACH2", "LDLR", "BMPR1B", "RFXAP", "ASPH", "PTK2B", "SLC1A5", "ENO2", "TRPM8", "SATB1", "MIER1", "SRSF1", "ATF3", "CCL5", "MCM6", "GCH1", "CAV1", "SLC20A1")
data=read.csv('Differentially Expressed mRNAs_down.csv')
entrezIDs=data$EntrezID
entrezIDs=entrezIDs[data$EntrezID %in% universe]
#
geneList=factor(as.integer(universe %in% entrezIDs))
names(geneList)=universe
GOdata=new('topGOdata',
           description='first try',
           ontology='CC',
           allGenes=geneList,
           nodeSize=5,
           annot=annFUN.org,
           mapping='org.Hs.eg.db',
           ID='entrez'
)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)
test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.stat)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)


allRes <- GenTable(GOdata, classic = resultFisher, KS = resultKS, weight = resultWeight,
                    orderBy = "weight", ranksOf = "classic", topNodes = 10)

#showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
gobarplot=function(allRes){
g=ggplot(allRes[1:10,],aes(y=-log10(as.numeric(weight)),x=Term))+
  geom_bar(stat='identity',fill='#56B4E9')+
  geom_text(aes(label=paste(Significant,'genes')),color='black',position=position_dodge(1),size=5,hjust=1.1)+
  theme_minimal()+coord_flip()+
  ylab('Enrichmen Score(-log10(pvalue))')+xlab('')+ggtitle('Down regulation CC')+
  theme(plot.title = element_text(hjust = 0.5,face='bold',size=16))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))
g
}
gobarplot(allRes)
