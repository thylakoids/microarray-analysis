rm(list=ls())
library("org.Hs.eg.db")
library("GSEABase")
library("GOstats")
library(RSQLite)
library(Category)
library("pathview")
enrich <- function(entrezIDs, orgDbName="org.Hs.eg.db", pvalueCutoff=.01){
  require(orgDbName, character.only=TRUE)
  require("GSEABase")
  require("GOstats")
  require("Category")
  goAnn <- get(gsub(".db", "GO", orgDbName))
  universe <- Lkeys(goAnn)
  onto <- c("BP", "MF", "CC")
  res <- lapply(onto, function(.onto){
    param <- new('GOHyperGParams',
                 geneIds= entrezIDs,
                 universeGeneIds=universe,
                 annotation=orgDbName,
                 ontology=.onto,
                 pvalueCutoff=pvalueCutoff,
                 conditional=FALSE,
                 testDirection="over")
    over <- hyperGTest(param)
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse=";")
    })
    summary <- summary(over)
    if(nrow(summary)>1) summary$Symbols <- glist[as.character(summary[, 1])]
    summary
  })
  names(res) <- onto
  keggAnn <- get(gsub(".db", "PATH", orgDbName))
  universe <- Lkeys(keggAnn)
  param <- new("KEGGHyperGParams",
               geneIds=entrezIDs,
               universeGeneIds=universe,
               annotation=orgDbName,
               categoryName="KEGG",
               pvalueCutoff=pvalueCutoff,
               testDirection="over")
  over <- hyperGTest(param)
  kegg <- summary(over)
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
  })
  kegg$Symbols <- glist[as.character(kegg$KEGGID)]
  res[["kegg"]] <- kegg
  res
}

genes <- c("AREG", "FKBP5", "CXCL13", "KLF9", "ZC3H12A", "P4HA1", "TLE1", "CREB3L2", "TXNIP", "PBX1", "GJA1", "ITGB8", "CCL3", "CCND2", "KCNJ15", "CFLAR", "CXCL10", "CYSLTR1", "IGFBP7", "RHOB", "MAP3K5", "CAV2", "CAPN2", "AKAP13", "RND3", "IL6ST", "RGS1", "IRF4", "G3BP1", "SEL1L", "VEGFA", "SMAD1", "CCND1", "CLEC3B", "NEB", "AMD1", "PDCD4", "SCD", "TM2D3", "BACH2", "LDLR", "BMPR1B", "RFXAP", "ASPH", "PTK2B", "SLC1A5", "ENO2", "TRPM8", "SATB1", "MIER1", "SRSF1", "ATF3", "CCL5", "MCM6", "GCH1", "CAV1", "SLC20A1")
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
go <- enrich(entrezIDs, "org.Hs.eg.db", .05)