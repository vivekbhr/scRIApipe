#!/usr/bin/env Rscript

## take the TCC matrix for each sample, subset it based on the merged EClist for that sample,
## then concatenate with the subsetted matrices from other samples
library(Matrix)
Args <- commandArgs(trailingOnly = TRUE)

matrixFileList <- Args[1] #c("transcripts_quant/ESC_merged/eq_counts/tcc.mtx", "transcripts_quant/NPC_merged/eq_counts/tcc.mtx")
ecMapList <- Args[2] #c("test/ESC_EC_subset_merged.txt", "test/NPC_EC_subset_merged.txt")
bclist <- Args[3] # "transcripts_quant/{sample}/eq_counts/tcc.barcodes.txt"
samples <- Args[4] # {sample}

out_mtx <- Args[5] # output merged mtx file
out_ec <- Args[6] # output EC-gene-txmap
out_bc <- Args[7] # output cell barcodes
matrixFileList <- trimws(unlist(strsplit(matrixFileList, ",", fixed = TRUE)))
ecMapList <- trimws(unlist(strsplit(ecMapList, ",", fixed = TRUE)))
bclist <- trimws(unlist(strsplit(bclist, ",", fixed = TRUE)))
samples <- trimws(unlist(strsplit(samples, ",", fixed = TRUE)))
## merge the EC map list based on the contained transcripts
eclist <- lapply(ecMapList, read.delim, header = FALSE,
                 col.names = c("EC", "Gene", "TxSet"), stringsAsFactors = FALSE)
eclist <- Reduce(function(x,y) merge(x, y, by = "TxSet"), eclist)

eclist1 <- eclist[ , grepl("EC", colnames(eclist))]
eclist2 <- cbind(eclist[,c("Gene.x", "TxSet")], eclist1)

## merge barcodes
bclist <- lapply(bclist, read.table, header = FALSE, stringsAsFactors = F)
bclist <- mapply(function(x,y) paste(x, y$V1, sep = "_"), samples, bclist)

## subset all matrices based on the merged EClist and merge them
subset_mtxList <- mapply(function(mtx,ec){
  ec_counts <- readMM(mtx)
  return(ec_counts[, ec + 1])# added 1 again, as the EC ids are 0-indexed
}, matrixFileList, eclist1)

writeMM(Reduce(rbind, subset_mtxList), file = out_mtx)
write.table(eclist2, file = out_ec, sep = "\t", row.names = F, col.names = F, quote = F)
write.table(unlist(bclist), file = out_bc, sep = "\n", row.names = F, col.names = F, quote = F)
