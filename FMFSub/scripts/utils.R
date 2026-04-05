suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

read_config <- function(path) yaml::read_yaml(path)

ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

read_tabular <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("tsv", "txt")) {
    fread(path, sep = "\t", data.table = FALSE)
  } else if (tolower(ext) %in% c("csv")) {
    fread(path, sep = ",", data.table = FALSE)
  } else if (tolower(ext) %in% c("xlsx", "xls")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Please install openxlsx to read Excel inputs.")
    }
    openxlsx::read.xlsx(path)
  } else {
    stop("Unsupported file extension: ", ext)
  }
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

coerce_group_factor <- function(meta, group_col, reference_group) {
  meta[[group_col]] <- factor(meta[[group_col]], levels = c(reference_group, setdiff(unique(meta[[group_col]]), reference_group)))
  meta
}

safe_model_terms <- function(meta, vars) {
  vars <- vars[vars %in% colnames(meta)]
  vars <- vars[!vapply(vars, function(v) all(is.na(meta[[v]])), logical(1))]
  vars
}

make_volcano <- function(res_df, lfc_col = "log2FoldChange", p_col = "padj", label_col = "gene",
                         lfc_thr = 1, alpha = 0.05, top_n = 15, title = NULL) {
  df <- res_df %>%
    mutate(.sig = !is.na(.data[[p_col]]) & .data[[p_col]] < alpha & abs(.data[[lfc_col]]) >= lfc_thr,
           .neglog10 = -log10(pmax(.data[[p_col]], 1e-300)))

  labels <- df %>%
    filter(.sig) %>%
    arrange(desc(.data[[lfc_col]])) %>%
    slice_head(n = top_n) %>%
    bind_rows(df %>% filter(.sig) %>% arrange(.data[[lfc_col]]) %>% slice_head(n = top_n)) %>%
    distinct(.data[[label_col]], .keep_all = TRUE)

  p <- ggplot(df, aes(x = .data[[lfc_col]], y = .neglog10)) +
    geom_point(aes(color = .sig), alpha = 0.7, size = 1.4) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
    geom_hline(yintercept = -log10(alpha), linetype = 2) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(x = "log2 fold-change", y = expression(-log[10](adjusted~p)), title = title)

  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(labels) > 0) {
    p <- p + ggrepel::geom_text_repel(data = labels, aes(label = .data[[label_col]]), size = 3)
  }
  p
}
