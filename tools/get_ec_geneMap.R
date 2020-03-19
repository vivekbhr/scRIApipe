#!/usr/bin/env Rscript

options("scipen"=100, "digits"=4)
### This script takes ~10 min  to process ~1.5 Mil lines
### get TCCs using bustools text output + create a EC to gene map

library(magrittr)
# test args
Args <- commandArgs(trailingOnly = TRUE)
tr2g <- Args[1]
# busfile <- Args[2]
transcript_list <- Args[2]
eclist <- Args[3]
ec_to_tr <- Args[4]
outFolder <- Args[5]
#nf <- as.numeric(Args[6]) # max number of fields in .ec file

# ## convert ec to tx map to ec to gene map
 tr2g <- read.delim(tr2g, header = FALSE, stringsAsFactors = FALSE)
 n <- c("tx_id", "gene_id")
 if(ncol(tr2g) == 4) {
   n <- c(n, "gene_name", "biotype")
 }
 colnames(tr2g) <- n

 tx_ids <- readr::read_delim(transcript_list, delim = "\t",
                             col_names = "tx", col_types = "c")
 tx_ids$gene_id <- tr2g[match(tx_ids$tx, tr2g$tx_id), "gene_id"]

getGenes <- function(x, pos, eclist = bus_eclist, file = out) {
  df <- dplyr::filter(x, ec %in% eclist)
  if(dim(df)[1] != 0) {
    apply(df, 1, function(x) {
      key <- as.integer(x[1])
      values <- unlist(strsplit(x[2], ","))
      values <- as.integer(values) + 1
      genes <- unique(tx_ids[values, ]$gene_id)
      name <- ifelse(length(genes) != 1, NA, genes)
      return(name)
    }) -> df$gene
    readr::write_tsv(df, file, append = TRUE)
    return(TRUE)
  } else {
    return(FALSE)
  }
}

bus_eclist <- readr::read_tsv(eclist, col_names="ec", col_types = "i")$ec
out <- file(file.path(outFolder, "ec-to-gene.txt"), open = "w")
readr::read_tsv_chunked(ec_to_tr, callback = readr::DataFrameCallback$new(getGenes), chunk_size = 50000,
                        col_names = c("ec", "tx"), col_types = "ic", trim_ws = TRUE, skip_empty_rows = TRUE)
close(out)
