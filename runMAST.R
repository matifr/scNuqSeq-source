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
  
  #Create data subsets for cortex and brainstem
  sca_cort <- subset(scaRaw, with(colData(scaRaw), B_region=='Cortex'))
  sca_bs <- subset(scaRaw, with(colData(scaRaw), B_region=='BrainStem'))
  
  print(assay(sca_cort))
  
  
  # calculate cellular detection rate
  cdr2 <-colSums(assay(sca_cort)>0)
  colData(sca_cort)$cngeneson <- scale(cdr2)
  cond<-factor(colData(sca_cort)$Condition)
  cond<-relevel(cond,"Control")
  colData(sca_cort)$Condition<-cond
  
  # calculate cellular detection rate
  cdr2_bs <-colSums(assay(sca_bs)>0)
  colData(sca_bs)$cngeneson <- scale(cdr2_bs)
  cond<-factor(colData(sca_bs)$Condition)
  cond<-relevel(cond,"Control")
  colData(sca_bs)$Condition<-cond
  
  print("Carrying out DE analysis..")
  zlm_sca_cort <- zlm(~Condition + cngeneson, sca_cort)
  
  #only test the cluster coefficient.
  summary_sca_cort <- summary(zlm_sca_cort, doLRT=TRUE)
  summary_Dt <- summary_sca_cort$datatable
  fcHurdle_cort <- merge(summary_Dt[contrast=='ConditionVE' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                         summary_Dt[contrast=='ConditionVE' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
  
  fcHurdle_cort <- fcHurdle_cort[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig_cort <- fcHurdle_cort[fdr<=0.05 ]
  setorder(fcHurdleSig_cort, fdr)
  
  
  # carry out DE analysis
  zlm_sca_bs <- zlm(~Condition + cngeneson, sca_bs)
  
  #only test the cluster coefficient.
  summary_sca_bs <- summary(zlm_sca_bs, doLRT=TRUE)
  summary_Dt <- summary_sca_bs$datatable
  fcHurdle_bs <- merge(summary_Dt[contrast=='ConditionVE' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                    summary_Dt[contrast=='ConditionVE' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
  
  fcHurdle_bs[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig_bs <- fcHurdle_bs[fdr<=0.05 ]
  print(dim(fcHurdleSig_bs))
  #fcHurdleSig_bs <- merge(fcHurdle_bs[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca_bs)), by='primerid')
  setorder(fcHurdleSig_bs, fdr)
  
  return( list("DEresults_cort" = fcHurdleSig_cort, "DEresults_bs" = fcHurdleSig_bs))
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


testMAST <- function(exp_filename, cdata_filename, save_filename_cort,save_filename_bs) {
  mylist <- loadData(exp_filename, cdata_filename)
  result <- runMAST(mylist)
  saveResult(result$DEresults_cort, save_filename_cort)
  saveResult(result$DEresults_bs, save_filename_bs)
}

# args should be:
# 1. exp_filename
# 2. cdata_filename
# 3. save_filename

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)

testMAST(args[1], args[2], args[3],args[4])





