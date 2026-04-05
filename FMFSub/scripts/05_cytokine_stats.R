suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
})
source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config <- read_config(args[1])
outdir <- file.path(config$output_dir, "cytokines")
ensure_dir(outdir)

cfg <- config$cytokines
x <- read_tabular(cfg$input_file)
meta <- read_tabular(cfg$metadata_file)

if (!all(c(cfg$sample_id_col, cfg$analyte_col, cfg$value_col) %in% colnames(x))) {
  sample_col <- intersect(c(cfg$sample_id_col, "sample_id"), colnames(x))[1]
  x <- x %>% pivot_longer(cols = -all_of(sample_col), names_to = cfg$analyte_col, values_to = cfg$value_col)
  names(x)[names(x) == sample_col] <- cfg$sample_id_col
}

dat <- x %>% left_join(meta, by = setNames(cfg$sample_id_col, cfg$sample_id_col))
if (isTRUE(cfg$log10_transform)) dat[[cfg$value_col]] <- log10(dat[[cfg$value_col]] + cfg$pseudocount)

stats <- dat %>%
  group_by(.data[[cfg$analyte_col]]) %>%
  wilcox_test(as.formula(paste(cfg$value_col, "~", cfg$group_col))) %>%
  adjust_pvalue(method = "BH") %>%
  arrange(p.adj)
write_tsv(stats, file.path(outdir, "cytokine_group_stats.tsv"))

sig_anal <- unique(stats[[cfg$analyte_col]])
pdf(file.path(outdir, "cytokine_boxplots.pdf"), width = 8, height = 6)
for (a in sig_anal) {
  dd <- dat %>% filter(.data[[cfg$analyte_col]] == a)
  print(ggplot(dd, aes_string(x = cfg$group_col, y = cfg$value_col, fill = cfg$group_col)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.6) +
          geom_jitter(width = 0.15, size = 1.5) + theme_bw(base_size = 12) +
          labs(title = a, x = NULL, y = ifelse(cfg$log10_transform, "log10 value", "value")) +
          theme(legend.position = "none"))
}
dev.off()
