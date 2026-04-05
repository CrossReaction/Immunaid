suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config <- read_config(args[1])
outdir <- file.path(config$output_dir, "somascan")
ensure_dir(outdir)

cfg <- config$somascan
expr <- read_tabular(cfg$expression_file)
meta <- read_tabular(cfg$metadata_file)
meta <- coerce_group_factor(meta, cfg$group_col, cfg$reference_group)
rownames(meta) <- meta[[cfg$sample_id_col]]

feature_col <- intersect(c(cfg$feature_id_col, "feature", "AptName", "Target"), colnames(expr))[1]
if (is.na(feature_col)) stop("Could not find feature identifier column in SomaScan matrix.")

sample_hits <- sum(colnames(expr) %in% rownames(meta))
if (sample_hits >= 2) {
  mat <- as.matrix(expr[, intersect(colnames(expr), rownames(meta)), drop = FALSE])
  rownames(mat) <- expr[[feature_col]]
} else {
  sample_col <- intersect(c(cfg$sample_id_col, "sample_id", "SampleID"), colnames(expr))[1]
  stopifnot(!is.na(sample_col))
  mat <- expr %>% column_to_rownames(sample_col) %>% as.matrix()
  mat <- t(mat)
}

common <- intersect(colnames(mat), rownames(meta))
mat <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
if (!isTRUE(cfg$log2_input)) mat <- log2(mat + 1)

covars <- safe_model_terms(meta, cfg$design_covariates)
design <- model.matrix(as.formula(paste("~ 0 +", paste(c(cfg$group_col, covars), collapse = " + "))), data = meta)
fit <- lmFit(mat, design)
contrast <- makeContrasts(contrasts = paste0(cfg$group_col, cfg$contrast_group, "-", cfg$group_col, cfg$reference_group), levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf, sort.by = "P") %>% rownames_to_column("feature")
write_tsv(res, file.path(outdir, "limma_results.tsv"))

pdf(file.path(outdir, "somascan_volcano.pdf"), width = 8, height = 6)
print(make_volcano(res %>% rename(log2FoldChange = logFC, padj = adj.P.Val, gene = feature),
                   lfc_thr = 0.4, alpha = cfg$alpha, top_n = 15,
                   title = paste(cfg$contrast_group, "vs", cfg$reference_group, "SomaScan")))
dev.off()
