#!/usr/bin/env Rscript

options("scipen"=100, "digits"=4)
####### Get Gene Counts per barcode by aggregating TCCs (and removing multi-genic TCCs)
library(Matrix)
library(magrittr)

#Args = c("test/ECtoGene_map.txt",
#         "test/output.mtx",
#         "test/output.ec.txt",
#         "test/output.barcodes.txt",
#         "test_gene")
Args <- commandArgs(trailingOnly = TRUE)
ec_to_gene <- Args[1]
mtx <- Args[2]
eclist <- Args[3]
bclist <- Args[4]
outFolder <- Args[5]

## load TCC output files
ec_genemap <- readr::read_tsv(ec_to_gene, col_names = c("ec", "gene"),
                              col_types = "ic")
bus_umi.mat <- Matrix::readMM(mtx)
barcodes <- readr::read_tsv(bclist, col_names = "bc", col_types = "c")$bc

## bus eclist and the bus_umi.mat still contains multi-genic ECs, but the ECtoGene map doesn't
## at this point it's best to discard these ECs
bus_eclist <- readr::read_tsv(eclist, col_names = "ec", col_types = "i")$ec
names(bus_eclist) <- ec_genemap[match(bus_eclist, ec_genemap$ec), ]$gene
rownames(bus_umi.mat) <- barcodes
colnames(bus_umi.mat) <- names(bus_eclist)
bus_umi.mat <- bus_umi.mat[,!is.na(colnames(bus_umi.mat))]
bus_genecounts <- fac2sparse(colnames(bus_umi.mat)) %*% t(bus_umi.mat) # returns gene * cell matrix

## write outputs : cell*gene matrix, barcodes, and genes
Matrix::writeMM(t(bus_genecounts), file = file.path(outFolder, "output.mtx"))
write.table(colnames(bus_genecounts), file = file.path(outFolder, "output.barcodes.txt"),
            row.names = F, col.names = F, quote = F)
write.table(rownames(bus_genecounts), file = file.path(outFolder, "output.genes.txt"),
            row.names = F, col.names = F, quote = F)
