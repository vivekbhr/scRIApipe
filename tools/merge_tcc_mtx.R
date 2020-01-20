#!/usr/bin/env Rscript

## take the TCC matrix for each sample, subset it based on the merged EClist for that sample,
## then concatenate with the subsetted matrices from other samples
library(Matrix)
Args <- commandArgs(trailingOnly = TRUE)

matrixFileList <- Args[1] #c("transcripts_quant/ESC_merged/eq_counts/output.mtx,transcripts_quant/NPC_merged/eq_counts/output.mtx")
ecMapList <- Args[2] #c("transcripts_quant/ESC_merged/eq_counts/ECtoGene_map.txt", "transcripts_quant/NPC_merged/eq_counts/ECtoGene_map.txt")
bclist <- Args[3] # "transcripts_quant/ESC_merged/eq_counts/output.barcodes.txt,transcripts_quant/NPC_merged/eq_counts/output.barcodes.txt"
samples <- Args[4] # "ESC,NPC"
mergeBy <- Args[5] # EC or TxSet

out_mtx <- Args[6] # output merged mtx file
out_ec <- Args[7] # output EC-gene-txmap
out_bc <- Args[8] # output cell barcodes


matrixFileList <- trimws(unlist(strsplit(matrixFileList, ",", fixed = TRUE)))
ecMapList <- trimws(unlist(strsplit(ecMapList, ",", fixed = TRUE)))
bclist <- trimws(unlist(strsplit(bclist, ",", fixed = TRUE)))
samples <- trimws(unlist(strsplit(samples, ",", fixed = TRUE)))

## merge barcodes
bclist <- lapply(bclist, read.table, header = FALSE, stringsAsFactors = F)
bclist <- lapply(seq_along(bclist), function(n) paste(samples[n], bclist[[n]]$V1, sep = "_"))


## merge the EC map list based on the contained transcripts
eclist <- lapply(ecMapList, read.delim, header = FALSE,
                 col.names = c("EC", "Gene", "TxSet"), stringsAsFactors = FALSE)

if(mergeBy == "EC") {
  eclist <- na.omit(Reduce(function(x,y) merge(x, y, by = "EC"), eclist))
  eclist2 <- eclist[1:3]
  subset_mtxList <- lapply(matrixFileList, function(mtx){
    ec_counts <- readMM(mtx)
    return(ec_counts[, eclist$EC + 1]) }) # added 1 again, as the EC ids are 0-indexed

} else {
  eclist <- na.omit(Reduce(function(x,y) merge(x, y, by = "TxSet"), eclist))
  eclist1 <- eclist[ , grepl("EC", colnames(eclist))]
  eclist2 <- cbind(eclist[,c("Gene.x", "TxSet")], eclist1)
  ## subset all matrices based on the merged EClist and merge them
  subset_mtxList <- mapply(function(mtx, ec){
    ec_counts <- readMM(mtx)
    return(ec_counts[, ec + 1])# added 1 again, as the EC ids are 0-indexed
  }, matrixFileList, eclist1)
}

## write outputs
writeMM(Reduce(rbind, subset_mtxList), file = out_mtx)
write.table(eclist2, file = out_ec, sep = "\t", row.names = F, col.names = F, quote = F)
write.table(unlist(bclist), file = out_bc, sep = "\n", row.names = F, col.names = F, quote = F)
