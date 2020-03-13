#!/usr/bin/env Rscript

library(Matrix)
library(ggplot2)
Args <- commandArgs(trailingOnly = TRUE)
## test Args
#Args <- c("../04b_scria_fullRun/transcripts_quant/TCCs_filtered_merged.mtx",
#          "../04b_scria_fullRun/transcripts_quant/barcodes_merged.txt",
#          "../04b_scria_fullRun/transcripts_quant/ECs_filtered_merged.txt",
#          "../04b_scria_fullRun/annotations/tr2g.tsv",
#          0.1, ".")

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
colnames(tcc.mtx) <- ecmap$V1
#cols <- colnames(tcc.mtx) == "ENSMUSG00000031575.18"
#colSums(tcc.mtx[,cols])

# keep track of associated txSet
kept_ecs <- colSums(tcc.mtx) > 10
tcc.mtx <- tcc.mtx[rowSums(tcc.mtx) > 50, kept_ecs]
#tcc.mtx <- tcc.mtx[rowSums(tcc.mtx) > 1,colSums(tcc.mtx) > 10]
ecmap <- ecmap[kept_ecs, ]

## logistic-reg + LRT per gene
lrt_gene <- function(gene, tcc_mtx, nullCompare = 2) {
  mtx <- as.data.frame(as.matrix(tcc_mtx[ , grepl(gene, colnames(tcc_mtx))])) # subset mtx by gene
  n <- paste("EC", 1:ncol(mtx), sep = "_") # col ID for the formula
  colnames(mtx) <- n
  mtx$label <- as.factor(gsub("(.*)_[AGTC]*", "\\1", rownames(mtx)))
  fmla <- as.formula(paste("label ~ ", paste(n, collapse= "+")))
  glm.fit <- glm(fmla, data = mtx, family = binomial)

  if(nullCompare == 1) {
    # this null is to get all genes
    fit_null <- glm("label ~ 1", data = mtx, family = binomial)
  } else {
    mtx$total <- rowSums(mtx[-ncol(mtx)])
    # this null is to compare to total count
    fit_null <- glm("label ~ total", data = mtx, family = binomial)
  }
  #summary(glm.fit)
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
