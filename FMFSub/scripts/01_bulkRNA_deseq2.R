suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(apeglm)
  library(pheatmap)
  library(EnhancedVolcano)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
})
source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "utils.R"))

args <- commandArgs(trailingOnly = TRUE)
config <- read_config(args[1])
outdir <- file.path(config$output_dir, "bulk_rnaseq")
ensure_dir(outdir)

cfg <- config$bulk_rnaseq
meta <- read_tabular(cfg$metadata_file)
meta <- coerce_group_factor(meta, cfg$group_col, cfg$reference_group)
rownames(meta) <- meta[[cfg$sample_id_col]]

if (file.exists(cfg$counts_file)) {
  counts_df <- read_tabular(cfg$counts_file)
  gene_col <- intersect(c("gene", "gene_id", "gene_symbol", "Gene"), colnames(counts_df))[1]
  stopifnot(!is.na(gene_col))
  count_mat <- as.matrix(counts_df[, setdiff(colnames(counts_df), gene_col), drop = FALSE])
  rownames(count_mat) <- counts_df[[gene_col]]
} else {
  tx2gene <- read_tabular(cfg$tx2gene_file)
  quant_dir <- cfg$quant_dir
  files <- file.path(quant_dir, meta[[cfg$sample_id_col]], "abundance.h5")
  if (!all(file.exists(files))) {
    files <- file.path(quant_dir, meta[[cfg$sample_id_col]], "abundance.tsv")
  }
  names(files) <- meta[[cfg$sample_id_col]]
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene[, 1:2])
  count_mat <- round(txi$counts)
}

common <- intersect(colnames(count_mat), rownames(meta))
count_mat <- count_mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

keep <- rowSums(count_mat >= cfg$min_count) >= ceiling(cfg$min_prop_samples * ncol(count_mat))
count_mat <- count_mat[keep, , drop = FALSE]

covars <- safe_model_terms(meta, cfg$design_covariates)
form <- as.formula(paste("~", paste(c(covars, cfg$group_col), collapse = " + ")))

dds <- DESeqDataSetFromMatrix(countData = round(count_mat), colData = meta, design = form)
dds <- DESeq(dds)
res <- results(dds, contrast = c(cfg$group_col, cfg$contrast_group, cfg$reference_group), alpha = cfg$alpha)
if (isTRUE(cfg$use_apeglm)) {
  coef_name <- grep(paste0("^", cfg$group_col, "_"), resultsNames(dds), value = TRUE)
  coef_name <- coef_name[grepl(cfg$contrast_group, coef_name) & grepl(cfg$reference_group, coef_name)]
  if (length(coef_name) == 1) {
    res <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
  }
}
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% arrange(padj)
write_tsv(res_df, file.path(outdir, "deseq2_results.tsv"))

sig_df <- res_df %>% filter(!is.na(padj), padj < cfg$alpha, abs(log2FoldChange) >= cfg$lfc_display_threshold)
write_tsv(sig_df, file.path(outdir, "deseq2_results_significant.tsv"))

vsd <- vst(dds, blind = FALSE)
write_tsv(as.data.frame(assay(vsd)) %>% rownames_to_column("gene"), file.path(outdir, "vst_matrix.tsv"))

pdf(file.path(outdir, "bulkRNA_qc_plots.pdf"), width = 8, height = 6)
plotPCA(vsd, intgroup = cfg$group_col) + ggtitle("VST PCA")
print(last_plot())
sample_dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dists), main = "Sample distance")
vol <- make_volcano(res_df, lfc_thr = cfg$lfc_display_threshold, alpha = cfg$alpha, top_n = cfg$top_label_n,
                    title = paste(cfg$contrast_group, "vs", cfg$reference_group))
print(vol)
dev.off()

message("Bulk RNA-seq analysis completed: ", outdir)
