#!/usr/bin/env Rscript

Args <- commandArgs(trailingOnly = TRUE)
genome <- Args[1]
gtf <- Args[2]
readLength <- as.numeric(Args[3])
outdir <- Args[4]

if (!require(BUSpaRse)) {
  devtools::install_github("BUStools/BUSpaRse", dep = FALSE)
}

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

g <- get(g)
GenomeInfoDb::seqlevelsStyle(g) <- "ensembl"

## txdb from gtf
message("Preparing annotations")
gtf.txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
GenomeInfoDb::seqlevelsStyle(gtf.txdb) <- "ensembl"
AnnotationDbi::saveDb(gtf.txdb, file.path(outdir, "gtf.txdb"))

## write cDNA fasta
message("Writing cDNA fasta")
txExons <- GenomicFeatures::exonsBy(gtf.txdb, by = "tx")
names(txExons) <- GenomicFeatures::transcripts(gtf.txdb)$tx_name
txSeqs <- GenomicFeatures::extractTranscriptSeqs(g, txExons)
Biostrings::writeXStringSet(txSeqs, filepath = file.path(outdir, "cDNA.fa"))

## create velocyto files (cDNA_introns.fasta, tr2g)
message("Writing velocity files")
BUSpaRse::get_velocity_files(gtf.txdb,
                             L = readLength - 1,
                             Genome = g,
                             out_path = outdir,
                             isoform_action = "collapse")
