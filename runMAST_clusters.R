#!/usr/bin/env Rscript
suppressMessages(library(MAST))
suppressPackageStartupMessages({library(data.table)})
suppressPackageStartupMessages({library(dplyr)})
options(mc.cores = 1) # gives me error messages when I use > 1

# Load the expression matrix and the cdata matrix
loadData <- function(exp_mat, cdata) {
  
  print("Loading data...")
  exp_mat <- (read.csv(exp_mat, row.names=NULL))
  rownames(exp_mat) <- exp_mat$cell_barcodes
  
  cdata <- (read.csv(cdata, row.names=NULL))
  rownames(cdata) <- cdata$cell_barcodes
  
  print("Finished loading data...")
  return(list('exp_mat' = exp_mat, 'cdata' =  cdata))
}


runMAST <- function(list) {
  
  print("Running MAST..")
  expr_matrix <- list$exp_mat
  expr_matrix <- expr_matrix %>% select(2:ncol(expr_matrix))
  
  cdata <- list$cdata
  fdata <- data.frame(primerid = names(expr_matrix), symbolid = names(expr_matrix))
  rownames(fdata) <- fdata$symbolid
  
  print(dim( t(expr_matrix)))
  print(dim(cdata))
  
  #print(head(t(expr_matrix)))
  #print(head(cdata))
  
  print("Make scaRaw Obj..")
  scaRaw <- FromMatrix(exprsArray = t(expr_matrix), cData = cdata, fData = fdata)
  
  # calculate cellular detection rate
  cdr2 <-colSums(assay(scaRaw)>0)
  colData(scaRaw)$cngeneson <- scale(cdr2)
  cond<-factor(colData(scaRaw)$groups)
  cond<-relevel(cond,"group_restCTRL")
  colData(scaRaw)$groups<-cond
  
  print("Carrying out DE analysis..")
  zlm_sca <- zlm(~groups + cngeneson, scaRaw)
  
  #only test the cluster coefficient.
  summary_sca <- summary(zlm_sca, doLRT=TRUE)
  summary_Dt <- summary_sca$datatable
  fcHurdle <- merge(summary_Dt[contrast=='groupsgroup_0_3_4VE' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                         summary_Dt[contrast=='groupsgroup_0_3_4VE' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
  
  fcHurdle <- fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig_cort <- fcHurdle[fdr<=0.05 ]
  setorder(fcHurdle, fdr)
  
  
  return( list("DEresults" = fcHurdle))
}

saveResult <- function(result, filename) {
  
  print("Saving results..")
  resultDf <- as.data.frame(result)
  colnames(resultDf)[1] = 'gene'
  colnames(resultDf)[2] = 'p'
  colnames(resultDf)[3] = 'logFC'
  colnames(resultDf)[6] = 'p.fdr.adj'
  resultDf <- resultDf[,c('gene','p','p.fdr.adj','logFC')]
  write.table(resultDf, file = filename, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
}


testMAST <- function(exp_filename, cdata_filename, save_filename,save_filename_bs) {
  mylist <- loadData(exp_filename, cdata_filename)
  result <- runMAST(mylist)
  saveResult(result$DEresults, save_filename)
}

# args should be:
# 1. exp_filename
# 2. cdata_filename
# 3. save_filename

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

testMAST(args[1], args[2], args[3])





