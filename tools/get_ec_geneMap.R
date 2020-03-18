#!/usr/bin/env Rscript

options("scipen"=100, "digits"=4)
### This script takes ~10 min  to process ~1.5 Mil lines
### get TCCs using bustools text output + create a EC to gene map

library(magrittr)
# test args
#Args = c("annotations/tr2g.tsv",
#         "transcripts_quant/ESC_merged/output.txt",
#         "transcripts_quant/ESC_merged/transcripts.txt",
#         "transcripts_quant/ESC_merged/matrix.ec",
#         "test")
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
#
#
# ## bus ops
# bus <- readr::read_delim(busfile, delim = "\t",
#                          col_names = c("bc", "umi", "ec", "count"),
#                          col_types = c("ffii")
#                          )
#
# ## prepare data frame of stats
# bus %>% dplyr::group_by(bc) %>% dplyr::summarise(
#                                      umi = length(unique(umi)),
#                                      reads = sum(count),
#                                      unique_reads = sum(count >= 1)) -> umi_perbc
# statsdf <- data.frame(quantile = c("Min", "1st_qt", "Median", "Mean", "3rd_qt", "Max"),
#                    nUMI = as.integer(summary(umi_perbc$umi)),
#                    nReads = as.integer(summary(umi_perbc$reads)),
#                    nDedup_Reads = as.integer(summary(umi_perbc$unique_reads))
#                    )
#
# ## unique reads are much more than nUMI..!!
# ## simple deupd per eq class
# split(bus, bus$bc) %>%
#   lapply(function(x) dplyr::group_by(x, ec) %>%
#            dplyr::summarise(umi_counts = sum(count >= 1))) %>%
#               plyr::ldply(data.frame) %>% tibble::as_tibble() -> bus_umi
# colnames(bus_umi)[1] <- "bc"
# bus_eclist <- sort(unique(bus_umi$ec))
#
# ## sanity check again
# # bus_umi %>% group_by(bc) %>% summarise(unique_reads = sum(umi_counts)) %>% summary(.$unique_reads)
# ## median 42K counts per cell!!
# #bus_umi_two <- tidyr::spread(bus_umi, bc, umi_counts, fill = 0)
# bus_umi %<>% tidyr::spread(ec, umi_counts, fill = 0)
# bus_umi.mat <- Matrix::Matrix(as.matrix(bus_umi[,2:ncol(bus_umi)]), sparse = TRUE)
#
# ## write output mtx, EC, barcodeList, stats
# Matrix::writeMM(bus_umi.mat, file = file.path(outFolder, "output.mtx"))
# write.table(bus_umi$bc, file = file.path(outFolder, "output.barcodes.txt"),
#             row.names = F, col.names = F, quote = F)
# write.table(bus_eclist, file = file.path(outFolder, "output.ec.txt"),
#             row.names = F, col.names = F, quote = F)
# write.table(statsdf, file = file.path(outFolder, "stats.txt"),
#             sep = "\t", quote = F, row.names = F)
#
#
# # convert eclist to char , to match with txid file
# #bus_eclist %<>% as.character()
#
# ## read an EC matrix and map EC to genes
# # if nf is given, skip the calculation (memory saving)
# #if(is.na(nf)) {
# #} else {
# #  no_col <- nf + 1
# #}
#
# ## map EC to genes,
# ## 1. remove ECs that don't map to genes (sanity check)
# ## 2. also remove ECs that map to multiple genes
# ## 3. also remove ECs that are not there in the bus text file
#
# #get_ec_gene_df <- function(line, outFile) {
# #  ecs <- read.table(text = gsub(",", "\t", line), fill = TRUE, header = F,
#                     #col.names = paste("Col", 1:no_col, sep = "_"),
# #                    sep = "\t", row.names = 1)
#   ## pre-filter ecs which are not in the ecList from the bus file
# #  ecs <- ecs[rownames(ecs) %in% bus_eclist, , drop = FALSE]
# #  if(length(ecs) > 0) {
#     # for each row in tx IDs
# #    apply(ecs, 1, function(x) {
# #      y <- as.numeric(na.omit(x)) + 1
#       # get corresponding geneID
# #      li <- na.omit(tx_ids[y, ]$gene_id)
#       # report NA if multiple genes are found
# #      li <- ifelse(length(unique(li)) > 1, NA, li)
# #      return(li)
# #    }) -> gene_ids
#
# #    which(!is.na(gene_ids)) -> kept_ec
# #    gene_ids[!is.na(gene_ids)] -> kept_genes
# #    cat(paste(as.numeric(names(kept_ec)), gene_ids[kept_ec], sep = "\t"), sep = "\n", append = TRUE, file = outFile)
# #  }
# #}
#
# #con <- file(ec_to_tr)
#
# #line <- c()
# #while(TRUE) {
# #  line = readLines(con, 50000)
# #  if(length(line) == 0) break
# #  else {
# #    cf <- count.fields(file(ec_to_tr), sep = ",")
# #    no_col <- max(cf) + 1
# #    get_ec_gene_df(line, out)
# #  }
# #}
#
# #df <- readr::read_tsv(ec_to_tr, col_names = c("ec", "tx"), col_types = "ic")

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
