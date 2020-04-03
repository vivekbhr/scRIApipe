#!/usr/bin/env Rscript
#setwd("/hpc/hub_oudenaarden/vbhardwaj/2019_vasa_seq/mESC_NPC/03_kallisto/05_scria_output")

library(Matrix)
library(ggplot2)

## test Args
#Args <- c("transcripts_quant/TCCs_filtered_merged.mtx",
#          "transcripts_quant/barcodes_merged.txt",
#          "transcripts_quant/ECs_filtered_merged.txt",
#          "annotations/tr2g.tsv",
#          0.1, "DTU_testing/ES_NPC_logisticReg")

## INPUT ARGS
Args <- commandArgs(trailingOnly = TRUE)

tcc.mtx <- readMM(Args[1])
bc <- read.table(Args[2], header = FALSE, stringsAsFactors = FALSE)$V1
ecmap <- read.delim(Args[3], header = FALSE, stringsAsFactors = FALSE)
tr2g <- read.table(Args[4], header = FALSE, stringsAsFactors = FALSE)
padj_threshold <- as.numeric(Args[5])

## OUTPUT ARGS
bcClusterMap <- read.delim(Args[6], header = FALSE, stringsAsFactors = FALSE)
threads <- Args[7]
outprefix <- Args[8]


## --------------- Functions ---------------------

## get DTU
lrt_gene <- function(gene, tcc_mtx, labels) {
  mtx <- as.data.frame(as.matrix(tcc_mtx[ , grepl(gene, colnames(tcc_mtx))])) # subset mtx by gene
  n <- paste("EC", 1:ncol(mtx), sep = "_") # col ID for the formula
  colnames(mtx) <- n
  mtx$label <- labels
  fmla <- as.formula(paste("label ~ ", paste(n, collapse= "+")))
  glm.fit <- glm(fmla, data = mtx, family = binomial)
  fit_null <- glm("label ~ 1", data = mtx, family = binomial)
  summary(glm.fit)
  #glm.probs <- predict(glm.fit,type = "response")
  lra <- anova(glm.fit, fit_null, test = "Chisq")
  return(lra$`Pr(>Chi)`[2])
}


## plot ECs of a given gene
plot_gene <- function(gene, ec_map, tcc_mtx, labels) {
  print(gene)
  test <- as.matrix(tcc_mtx[ , grepl(gene, colnames(tcc_mtx)), drop = FALSE])
  ecmap2 <- ec_map[grepl(gene, ec_map$V1), ] # subset ecmap by gene
  n <- paste0("TxSet: ", ecmap2$V2) # use the tx set as col IDs
  #colnames(test) <- paste0("EC_", 1:ncol(test))
  colnames(test) <- n
  t <- as.data.frame(test)
  t$labels <- labels
  t <- reshape2::melt(t)
  ggplot(t, aes(labels, value, col = labels)) + geom_jitter(height = 0.01) +
    labs(x = "Sample", y = "log10(Counts)", title = gene) + facet_wrap(~variable) +
    scale_y_log10()
}


## ---------------

## --------------- Prepare --------------
# rownames = sample_barcodes (unique)
rownames(tcc.mtx) <- bc
# colnames = gene name (non-unique)
colnames(tcc.mtx) <- ecmap$V1
#cols <- colnames(tcc.mtx) == "ENSMUSG00000031575.18"

# Filter TCC, but keep track of associated txSet/labels
kept_ecs <- colSums(tcc.mtx) > 10
kept_bcs <- rowSums(tcc.mtx) > 50

tcc.mtx <- tcc.mtx[kept_bcs, kept_ecs]
ecmap <- ecmap[kept_ecs, ]
bcLabels <- bcClusterMap[kept_bcs, ]


## ------------- Execute -------------

genes <- unique(colnames(tcc.mtx))
## parallel
system.time(
  plist <- unlist(mclapply(genes, function(x) lrt_gene(x, tcc.mtx, bcLabels), 
                           mc.cores = threads))
)

if(length(plist) == 0) {
  warning("Output empty!!")
  quit(save = "no", status = 1, runLast = FALSE)
} else {
  out <- data.frame(gene = genes, pval = plist)
}

out$padj <- p.adjust(out$pval, method = "BH")
sigGenes <- as.character(out[out$padj < padj_threshold, ]$genes)
topGenes <- as.character(out[order(out$padj, decreasing=F), "genes"][1:20])
message(paste0(length(sigGenes), " significant genes left after FDR correction!!"))


## ------------- Save Output -------------

pdf(paste0(outprefix, "_sigGenes_plots.pdf"))
lapply(topGenes, plot_gene, ec_map = ecmap)
dev.off()
## get names for sigGenes
gn <- unique(tr2g[match(sigGenes, tr2g$V2), "V3"])
write.table(gn, file = paste0(outprefix, "_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigGenes, file = paste0(outprefix, "_sigGenes_ensID.txt"), quote = F, row.names = F, col.names = F)


############### --------------- UNUSED --------------------  ###########
### pca/umap
#tcc.mtx2 <- log2((tcc.mtx/colSums(tcc.mtx))*100000 + 0.1)
#tcc.mtx2 <- tfidf(t(tcc.mtx))
#tcc.pc <- irlba::irlba(tcc.mtx2, nv = 50)
#pc <- tcc.pc$u %*% diag(tcc.pc$d) %>% as.data.frame()
#colnames(pc) <- paste0("PC", 1:ncol(pc))
#pc <- do_pca(tcc.mtx2) %>% as.data.frame()
#n <- as.factor(gsub("(.*)_[AGTC]*", "\\1", rownames(tcc.mtx)))
#pc$cell <- n
#ggplot(pc, aes(PC1, PC2, col = cell)) + geom_point()

#uwot::umap(pc[1:50], spread = 5, n_neighbors = 15, min_dist = 0.1) %>% as.data.frame() -> final.umap
#colnames(final.umap) <- c("UMAP1", "UMAP2")
#final.umap$cell <- n
#final.umap$cell <- as.factor(gsub("(ESC|NPC)_.*", "\\1", rownames(tcc.mtx)))
#ggplot(final.umap, aes(UMAP1, UMAP2, col = cell)) + geom_point()

##
#tr <- read.delim("transcripts_quant/ESC_merged/transcripts.txt", header = F)
## Second approach :: DEXSeq followed by lancaster p-value aggregation
#library(DEXSeq)
#exnames <- colnames(tcc.mtx)

#sapply(unique(exnames), function(x){
#  erep <- sum(grepl(x, exnames, fixed = TRUE))
#  erep <- paste0("E", 1:erep)
#  return(erep)
#}) %>% unlist() -> featureID

#cd <- t(tcc.mtx)
#rownames(cd) <- exnames2

#dds <- DEXSeqDataSet(countData = as.matrix(cd),
#              sampleData = data.frame(row.names = rownames(tcc.mtx),
#                                      condition = gsub("(.*)_[AGTC]*", "\\1", rownames(tcc.mtx))
#                                      ),
#              groupID = exnames,
#              featureID = featureID)
#dds %<>% estimateSizeFactors()
