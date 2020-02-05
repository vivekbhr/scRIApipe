#!/usr/bin/env Rscript
#setwd("/hpc/hub_oudenaarden/vbhardwaj/2019_vasa_seq/mESC_NPC/03_kallisto/05_scria_output")

library(Matrix)
library(ggplot2)
Args <- commandArgs(trailingOnly = TRUE)
## test Args
#Args <- c("transcripts_quant/TCCs_filtered_merged.mtx",
#          "transcripts_quant/barcodes_merged.txt",
#          "transcripts_quant/ECs_filtered_merged.txt",
#          "annotations/tr2g.tsv",
#          0.1, "DTU_testing/ES_NPC_logisticReg")

## INPUT ARGS
tcc.mtx <- readMM(Args[1])
bc <- read.table(Args[2], header = FALSE, stringsAsFactors = FALSE)$V1
ecmap <- read.delim(Args[3], header = FALSE, stringsAsFactors = FALSE)
tr2g <- read.table(Args[4], header = FALSE, stringsAsFactors = FALSE)
padj_threshold <- as.numeric(Args[5])

## OUTPUT ARGS
outprefix <- Args[6]

#labels <- gsub("(.*)_[AGTC]*", "\\1", bc) ## need to find a way to use externally defined label
# rownames = sample_barcodes (unique)
rownames(tcc.mtx) <- bc
# colnames = gene name (non-unique)
colnames(tcc.mtx) <- ecmap$V2
#cols <- colnames(tcc.mtx) == "ENSMUSG00000031575.18"
#colSums(tcc.mtx[,cols])

# keep track of associated txSet
kept_ecs <- colSums(tcc.mtx) > 10
tcc.mtx <- tcc.mtx[rowSums(tcc.mtx) > 50, kept_ecs]
#tcc.mtx <- tcc.mtx[rowSums(tcc.mtx) > 1,colSums(tcc.mtx) > 10]
ecmap <- ecmap[kept_ecs, ]


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

## logistic-reg + LRT per gene
lrt_gene <- function(gene, tcc_mtx) {
  mtx <- as.data.frame(as.matrix(tcc_mtx[ , grepl(gene, colnames(tcc_mtx))])) # subset mtx by gene
  n <- paste("EC", 1:ncol(mtx), sep = "_") # col ID for the formula
  colnames(mtx) <- n
  mtx$label <- as.factor(gsub("(.*)_[AGTC]*", "\\1", rownames(mtx)))
  fmla <- as.formula(paste("label ~ ", paste(n, collapse= "+")))
  glm.fit <- glm(fmla, data = mtx, family = binomial)
  fit_null <- glm("label ~ 1", data = mtx, family = binomial)
  summary(glm.fit)
  #glm.probs <- predict(glm.fit,type = "response")
  lra <- anova(glm.fit, fit_null, test = "Chisq")
  return(lra$`Pr(>Chi)`[2])
}

genes <- unique(colnames(tcc.mtx))
out <- data.frame(genes = genes,
                  pval = sapply(genes, lrt_gene, tcc_mtx = tcc.mtx))

if(dim(out)[1] == 0) {
  warning("Output empty!!")
}

out$padj <- p.adjust(out$pval, method = "BH")
sigGenes <- as.character(out[out$padj < padj_threshold, ]$genes)

message(paste0(length(sigGenes), " significant genes left after FDR correction!!"))

## plot a given gene
plot_gene <- function(gene, ec_map) {
  print(gene)
  test <- as.matrix(tcc.mtx[ , grepl(gene, colnames(tcc.mtx)), drop = FALSE])
  ecmap2 <- ec_map[grepl(gene, ec_map$V2), ] # subset ecmap by gene
  n <- paste0("TxSet: ", ecmap2$V3) # use the tx set as col IDs
  #colnames(test) <- paste0("EC_", 1:ncol(test))
  colnames(test) <- n
  t <- as.data.frame(test)
  t$labels <- gsub("(.*)_[AGTC]*", "\\1", rownames(t))
  t <- reshape2::melt(t)
  ggplot(t, aes(labels, value, col = labels)) + geom_jitter(height = 0.01) +
    labs(x = "Sample", y = "log10(Counts)", title = gene) + facet_wrap(~variable) +
    scale_y_log10()
}

pdf(paste0(outprefix, "_sigGenes_plots.pdf"))
lapply(sigGenes, plot_gene, ec_map = ecmap)
dev.off()

## get names for sigGenes
gn <- unique(tr2g[match(sigGenes, tr2g$V2), "V3"])
write.table(gn, file = paste0(outprefix, "_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigGenes, file = paste0(outprefix, "_sigGenes_ensID.txt"), quote = F, row.names = F, col.names = F)

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
