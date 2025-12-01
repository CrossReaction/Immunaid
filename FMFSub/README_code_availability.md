# Code Availability – RP / FMF Multi-omics

This repository contains reproducible analysis code for three data types used in the manuscript:

1. **Bulk RNA-seq** — differential expression, normalization and downstream reporting.
2. **SomaScan proteomics** — differential abundance analysis using limma.
3. **Luminex/ELISA cytokines** — non-parametric testing with multiple-testing correction, effect sizes, and visualization.

## Contents

- `gene_cpm_deseq2_fmfnc_20250913_EN.ipynb` — Python notebook for bulk RNA-seq CPM-based workflows (pyDESeq2).
- `SomaScan_disease_nc_limma_analysis_fmfnc_20250914_EN.ipynb` — R-based notebook for SomaScan analysis with limma.
- `IPA_data_preparation_20250915_EN.ipynb` — Python notebook to assemble multi-omics tables for IPA and perform cytokine stats.

> _EN_ files contain English-only annotations and comments suitable for public release.

## Environments

- **Python** ≥ 3.10  
  - `pandas`, `numpy`, `matplotlib`, `seaborn`, `statsmodels`, `scipy`, `pyDESeq2`, `gprofiler-official`

- **R** ≥ 4.2  
  - `limma`, `edgeR`, `tidyverse`, `readr`, `stringr`

> See `environment.yml` and `renv.lock` (optional) for exact versions.

## Reproducibility

- All analysis is performed at the **sample level** (one row = one biological replicate).
- Differential analyses include FDR control (Benjamini–Hochberg).
- Figures are produced with fixed random seeds where applicable.

## How to run

```bash
# option 1: open notebooks
jupyter lab

# option 2: convert notebooks to scripts
jupyter nbconvert --to script *.ipynb
```

## License

MIT © 2025 Zhicheng Zhou@Crossreaction
