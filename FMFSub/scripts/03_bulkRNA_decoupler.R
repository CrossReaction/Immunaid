suppressPackageStartupMessages({
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config <- read_config(args[1])
outdir <- file.path(config$output_dir, "bulk_rnaseq", "decoupler")
ensure_dir(outdir)

vst_df <- read_tabular(file.path(config$output_dir, "bulk_rnaseq", "vst_matrix.tsv"))
meta <- read_tabular(config$bulk_rnaseq$metadata_file)
rownames(meta) <- meta[[config$bulk_rnaseq$sample_id_col]]
mat <- vst_df %>% column_to_rownames("gene") %>% as.matrix()
mat <- scale(mat)

if (!requireNamespace("OmnipathR", quietly = TRUE)) {
  stop("Please install OmnipathR for decoupleR resources.")
}
prog <- decoupleR::get_progeny(organism = 'human', top = 100)
doro <- decoupleR::get_dorothea(organism = 'human', levels = c('A','B'))

path_scores <- run_mlm(mat = mat, network = prog, .source = 'source', .target = 'target', .mor = 'weight')
tf_scores <- run_mlm(mat = mat, network = doro, .source = 'source', .target = 'target', .mor = 'mor')

write_tsv(path_scores, file.path(outdir, "progeny_scores_long.tsv"))
write_tsv(tf_scores, file.path(outdir, "dorothea_scores_long.tsv"))

plot_scores <- function(df, outfile, title) {
  df2 <- df %>%
    left_join(meta %>% rownames_to_column("condition"), by = "condition") %>%
    group_by(source, .data[[config$bulk_rnaseq$group_col]]) %>%
    summarize(score = mean(score, na.rm = TRUE), .groups = "drop")
  pdf(outfile, width = 10, height = 6)
  print(ggplot(df2, aes(x = reorder(source, score), y = score, fill = .data[[config$bulk_rnaseq$group_col]])) +
          geom_col(position = "dodge") + coord_flip() + theme_bw(base_size = 12) + labs(x = NULL, y = "Mean MLM score", title = title))
  dev.off()
}
plot_scores(path_scores, file.path(outdir, "progeny_barplot.pdf"), "PROGENy pathway activity")
plot_scores(tf_scores, file.path(outdir, "dorothea_barplot.pdf"), "DoRothEA TF activity")
