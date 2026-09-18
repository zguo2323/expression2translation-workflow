suppressPackageStartupMessages({library(tximport); library(DESeq2); library(jsonlite)})
task <- fromJSON(commandArgs(trailingOnly=TRUE)[1])
cfg <- task$config
out <- task$directory
dir.create(out, recursive=TRUE, showWarnings=FALSE)
samples <- task$samples
stopifnot(!anyDuplicated(samples$sample_id), all(file.exists(samples$quant)))
tx2gene <- read.delim(task$tx2gene, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(identical(names(tx2gene), c("transcript_id", "gene_id")), !anyDuplicated(tx2gene$transcript_id))
for (file in samples$quant) {
  q <- read.delim(file, check.names=FALSE)
  stopifnot(!anyDuplicated(q$Name), setequal(q$Name, tx2gene$transcript_id),
            all(is.finite(q$NumReads)), all(q$NumReads >= 0), all(q$EffectiveLength > 0))
}
files <- setNames(samples$quant, samples$sample_id)
# Estimated fragment counts plus sample-specific average transcript lengths.
# Never feed TPM or independently rounded abundance values to DESeq2.
txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no", dropInfReps=TRUE)
write_matrix <- function(matrix, filename, id_column="gene_id") {
  tab <- data.frame(id=rownames(matrix), matrix, check.names=FALSE)
  names(tab)[1] <- id_column
  write.table(tab, file=file.path(out, filename),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
write_matrix(txi$counts, "gene_estimated_counts.tsv")
write_matrix(txi$abundance, "gene_tpm.tsv")
write_matrix(txi$length, "gene_average_length.tsv")
saveRDS(txi, file.path(out, "tximport.rds"))
write.table(samples, file.path(out, "samples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
summary <- list(synthetic=cfg$synthetic, reference_id=task$reference_id, config=cfg,
                design="~ condition", countsFromAbundance="no", genes_imported=nrow(txi$counts),
                deseq2_run=FALSE, contrast=cfg$contrast,
                versions=list(R=R.version.string, tximport=as.character(packageVersion("tximport")),
                              DESeq2=as.character(packageVersion("DESeq2"))))
if (isTRUE(cfg$run_deseq2)) {
  stopifnot(isTRUE(cfg$design_confirmed), nzchar(cfg$design_evidence), all(table(samples$condition) >= 2))
  coldata <- data.frame(condition=factor(samples$condition, levels=c(cfg$contrast[3], cfg$contrast[2])), row.names=samples$sample_id)
  stopifnot(!anyNA(coldata$condition))
  dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)
  keep <- rowSums(counts(dds) >= cfg$min_count) >= cfg$min_samples
  stopifnot(sum(keep) >= 2)
  dds <- dds[keep, ]
  dds <- DESeq(dds, fitType=cfg$fit_type, parallel=FALSE, minReplicatesForReplace=Inf)
  res <- results(dds, contrast=cfg$contrast, alpha=cfg$alpha, independentFiltering=TRUE)
  write_matrix(counts(dds, normalized=TRUE), "gene_normalized_counts.tsv")
  table <- data.frame(gene_id=rownames(res), as.data.frame(res), check.names=FALSE)
  write.table(table, file.path(out, "differential_expression.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  # varianceStabilizingTransformation also works with fewer than 1000 genes.
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  expression <- assay(vsd)
  variable <- apply(expression, 1, var) > 0
  stopifnot(sum(variable) >= 2)
  pc <- prcomp(t(expression[variable, , drop=FALSE]), scale.=FALSE)
  pca <- data.frame(sample_id=rownames(pc$x), condition=coldata[rownames(pc$x), "condition"], pc$x[, 1:2, drop=FALSE])
  write.table(pca, file.path(out, "pca.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write_matrix(cor(expression), "sample_correlation.tsv", "sample_id")
  pdf(file.path(out, "pca.pdf"), width=6, height=5)
  plot(pca$PC1, pca$PC2, col=as.integer(pca$condition), pch=19, xlab="PC1", ylab="PC2",
       main=if (cfg$synthetic) "SYNTHETIC TEST: RNA PCA" else "RNA PCA (VST, blind=FALSE)")
  text(pca$PC1, pca$PC2, labels=pca$sample_id, pos=3, cex=0.65)
  dev.off()
  saveRDS(dds, file.path(out, "deseq2.rds"))
  summary$deseq2_run <- TRUE
  summary$genes_tested <- nrow(dds)
  summary$significant_genes <- sum(res$padj < cfg$alpha, na.rm=TRUE)
  summary$actual_fit_type <- attr(dispersionFunction(dds), "fitType")
  summary$positive_log2FC <- paste(cfg$contrast[2], "higher than", cfg$contrast[3])
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
summary$inputs_md5 <- as.list(tools::md5sum(c(samples$quant, task$tx2gene, "src/rnaseq/gene_analysis.R")))
write_json(summary, file.path(out, "analysis.json"), pretty=TRUE, auto_unbox=TRUE, na="null")
