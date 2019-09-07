#!/usr/bin/env Rscript

Args <- commandArgs(trailingOnly = TRUE)
genome <- Args[1]
gtf <- Args[2] #"../01_annotations/ens95_gene_annotations.gtf"
readLength <- as.numeric(Args[3])
outdir <- Args[4] #"../01_annotations/ens95_velocyto_collapse"

## install required packages
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(BiocManager)) {
  install.packages("BiocManager")
}
# Install from GitHub
#BiocManager::install(c("BSgenome", "GenomicFeatures", "AnnotationDbi", "AnnotationHub", "pcaMethods"))

if (!require(BUSpaRse)) devtools::install_github("BUStools/BUSpaRse")
if (!require(seurat-wrappers)) devtools::install_github("satijalab/seurat-wrappers")
if (!require(seurat-wrappers)) devtools::install_github("velocyto-team/velocyto.R")


## load appropriate bsgenome
message("Loading genome")
avail_genomes <- BSgenome::available.genomes()
g <- grep(genome, avail_genomes, value = TRUE)
g <- grep("masked", g, invert = TRUE, value = TRUE)

g_load <- suppressWarnings(require(g, character.only = TRUE))

if(!(isTRUE(g_load))) {
  BiocManager::install(g, update = FALSE)
  require(g, character.only = TRUE)
  }
#library(g)
g <- get(g)
GenomeInfoDb::seqlevelsStyle(g) <- "ensembl"

## txdb from gtf
message("Preparing annotations")
gtf.txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
AnnotationDbi::saveDb(gtf.txdb, file.path(outdir, "gtf.txdb"))

## create velocyto files
message("Writing velocity files")
BUSpaRse::get_velocity_files(gtf.txdb,
                   L = readLength - 1,
                   Genome = g,
                   out_path = outdir,
                   isoform_action = "collapse")
