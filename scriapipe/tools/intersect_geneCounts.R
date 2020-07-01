#!/usr/bin/env Rscript

####### Match spliced/unspliced (or Any 2) matrices by gene names (column), and write back only common entries
library(Matrix)
Args <- commandArgs(trailingOnly = TRUE)
splicedMtx <- Args[1]
unsplicedMtx <- Args[2]
splicedGenes <- Args[3]
unsplicedGenes <- Args[4]
splicedOut <- Args[5]
unsplicedOut <- Args[6]


## load  output files
mt <- Matrix::readMM(splicedMtx)
umt <- Matrix::readMM(unsplicedMtx)
colnames(mt) <- read.delim(splicedGenes, header = FALSE)$V1
colnames(umt) <- read.delim(unsplicedGenes, header = FALSE)$V1

is <- intersect(colnames(mt), colnames(umt))

## write outputs : cell*gene matrix, and (intersected) genes
Matrix::writeMM(mt[ ,is], file = file.path(splicedOut, "output_isect.mtx"))
Matrix::writeMM(umt[ ,is], file = file.path(unsplicedOut, "output_isect.mtx"))
write.table(is, col.names = F, row.names = F, quote = F, file = file.path(splicedOut, "output.genes_isect.txt"))
write.table(is, col.names = F, row.names = F, quote = F, file = file.path(unsplicedOut, "output.genes_isect.txt"))