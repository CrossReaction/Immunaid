suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})
source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config <- read_config(args[1])
outdir <- file.path(config$output_dir, "bulk_rnaseq", "go")
ensure_dir(outdir)
res <- read_tabular(file.path(config$output_dir, "bulk_rnaseq", "deseq2_results_significant.tsv"))
all_genes <- read_tabular(file.path(config$output_dir, "bulk_rnaseq", "deseq2_results.tsv"))$gene

sig_genes <- unique(res$gene)
if (length(sig_genes) < 5) stop("Too few significant genes for GO analysis.")

ego <- enrichGO(gene = sig_genes,
                universe = all_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_df <- as.data.frame(ego)
write_tsv(ego_df, file.path(outdir, "go_enrichment.tsv"))

pdf(file.path(outdir, "go_dotplot.pdf"), width = 10, height = 8)
print(dotplot(ego, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free_y"))
dev.off()
