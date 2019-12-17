#!/usr/bin/env Rscript

options("scipen"=100, "digits"=4)
### This script takes ~10 min  to process ~1.5 Mil lines
Args <- commandArgs(trailingOnly = TRUE)

tr2g <- Args[1] # "annotations/tr2g.tsv"
transcript_list <- Args[2] # "transcripts_quant/{sample}/transcripts.txt"
ec_to_tr <- Args[3] # "transcripts_quant/{sample}/eq_counts/tcc.ec.txt"
tcc_mtx <- Args[4] # "transcripts_quant/{sample}/eq_counts/tcc.mtx"
out_ec_to_gene <- Args[5] # "EC to gene map"

## convert ec to tx map to ec to gene map
tr2g <- read.delim(tr2g, header = FALSE, stringsAsFactors = FALSE)
colnames(tr2g) <- c("tx_id", "gene_id", "gene_name", "biotype")

tx_ids <- read.delim(transcript_list, header = FALSE, stringsAsFactors = FALSE)
tx_ids$gene_id <- tr2g[match(tx_ids$V1, tr2g$tx_id), "gene_id"]

## read an EC matrix and map EC to genes
cf <- count.fields(file(ec_to_tr), sep = ",")
no_col <- max(cf) #- 1
#no_row <- length(cf)

## map EC to genes,
## 1. remove ECs that don't map to genes
## 2. also remove ECs that map to multiple genes
get_ec_gene_df <- function(line, outFile) {
  ecs <- read.table(text = gsub(",", "\t", line), fill = TRUE, header = F,
                    col.names = paste("Col", 1:no_col, sep = "_"), sep = "\t", row.names = 1)
  apply(ecs, 1, function(x) {
    y <- as.numeric(na.omit(x)) + 1
    li <- na.omit(tx_ids[y, ]$gene_id)
    li <- ifelse(length(unique(li)) > 1, NA, li)
    return(li)
  }) -> gene_ids

  which(!is.na(gene_ids)) -> kept_ec
  gene_ids[!is.na(gene_ids)] -> kept_genes
  cat(paste(as.numeric(names(kept_ec)), gene_ids[kept_ec], sep = "\t"), sep = "\n", append = TRUE, file = outFile)
}

con <- file(ec_to_tr, "r")
out <- file(out_ec_to_gene, open = "w")

lines <- c()
while(TRUE) {
  line = readLines(con, 50000)
  if(length(line) == 0) break
  else {
    get_ec_gene_df(line, out)
  }
}

close(out)
close(con)
