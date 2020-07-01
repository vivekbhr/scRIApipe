#!/usr/bin/env Rscript

## take the TCC matrix for each sample, subset it based on the merged EClist for that sample,
## then concatenate with the subsetted matrices from other samples
library(Matrix)
Args <- commandArgs(trailingOnly = TRUE)

matrixFileList <- Args[1] #"c("transcripts_quant/ESC_merged/eq_counts/output.mtx,transcripts_quant/NPC_merged/eq_counts/output.mtx")
ecList.filtered <- Args[2] #"transcripts_quant/ESC_merged/eq_counts/ECtoGene_map.txt", "transcripts_quant/NPC_merged/eq_counts/ECtoGene_map.txt"
ecList.unfiltered <- Args[3] #"transcripts_quant/ESC_merged/eq_counts/output.ec.txt", "transcripts_quant/NPC_merged/eq_counts/output.ec.txt"
bclist <- Args[4] # "transcripts_quant/ESC_merged/eq_counts/output.barcodes.txt,transcripts_quant/NPC_merged/eq_counts/output.barcodes.txt"
samples <- Args[5] # "ESC,NPC"
mergeBy <- Args[6] # EC or TxSet

out_mtx <- Args[7] # output merged mtx file
out_ec <- Args[8] # output EC-gene-txmap
out_bc <- Args[9] # output cell barcodes


matrixFileList <- trimws(unlist(strsplit(matrixFileList, ",", fixed = TRUE)))
ecList.filtered <- trimws(unlist(strsplit(ecList.filtered, ",", fixed = TRUE)))
ecList.unfiltered <- trimws(unlist(strsplit(ecList.unfiltered, ",", fixed = TRUE)))
bclist <- trimws(unlist(strsplit(bclist, ",", fixed = TRUE)))
samples <- trimws(unlist(strsplit(samples, ",", fixed = TRUE)))

## ---- if only one sample, just link the files, else merge them --- ##

## -----------merge barcodes
bclist <- lapply(bclist, read.table, header = FALSE, stringsAsFactors = F)
bclist <- lapply(seq_along(bclist), function(n) paste(samples[n], bclist[[n]]$V1, sep = "_"))
write.table(unlist(bclist), file = out_bc, sep = "\n", row.names = F, col.names = F, quote = F)


## ------------  merge EC and mtx list
if(length(matrixFileList) == 1 & length(ecList.filtered) == 1) {
  file.link(matrixFileList, out_mtx)
  ## reverse the colums (gene, tx, ec)
  eclist <- read.delim(ecList.filtered, header = FALSE,
            col.names = c("EC", "TxSet", "Gene"), stringsAsFactors = FALSE)
  write.table(eclist[c("Gene", "TxSet", "EC")], file = out_ec, sep = "\t", row.names = F, col.names = F, quote = F)
} else {
  
  ## merge the EC map list based on the contained transcripts or ECs (depending on "mergeBy")
  eclist <- lapply(ecList.filtered, read.delim, header = FALSE,
                   col.names = c("EC", "TxSet", "Gene"), stringsAsFactors = FALSE)
  
  ## first subset the mtx for each sample based on it's own filtered ECs
  mtxlist <- lapply(seq_along(matrixFileList), function(n){
    ec_counts <- readMM(matrixFileList[n])
    ecIDs <- read.delim(ecList.unfiltered[n], header = F)$V1
    keptEC <- ecIDs %in% eclist[[n]]$EC
    outmtx <- ec_counts[, keptEC]
    colnames(outmtx) <- ecIDs[keptEC]
    return(outmtx)
  })
  
  
  ## then subset for ECs that match across samples (dropping ECs specific to only one sample)
  if(mergeBy == "EC") {
    eclist.match <- na.omit(plyr::join_all(eclist, by = "EC", type = "inner"))
    eclist2 <- eclist.match[1:3]
    subset_mtxList <- lapply(seq_along(mtxlist), function(n){
      mt <- mtxlist[[n]]
      keptEC <- colnames(mt) %in% as.character(eclist.match$EC)
      return(mt[, keptEC]) })
    
  } else {
    eclist.match <- na.omit(plyr::join_all(eclist, by = "TxSet", type = "inner"))
    eclist1 <- eclist.match[ , grepl("EC", colnames(eclist.match))]
    eclist2 <- cbind(eclist.match[,c("Gene", "TxSet")], eclist1)
    ## subset all matrices based on the merged EClist and merge them (keeps more ECs than the one above)
    subset_mtxList <- lapply(seq_along(mtxlist), function(n){
      mt <- mtxlist[[n]]
      keptEC <- colnames(mt) %in% as.character(eclist1[, n])
      return(mt[, keptEC])
    })
  }
  
  ## write outputs
  writeMM(Reduce(rbind, subset_mtxList), file = out_mtx)
  write.table(eclist2, file = out_ec, sep = "\t", row.names = F, col.names = F, quote = F)
  
}


