# FMFSub multi-omics reproducibility code

Reproducible analysis code for the ImmunAID FMF genotype-dosage study comparing **biallelic FMF** (`n = 27`) versus **simple heterozygous FMF** (`n = 11`) during acute crisis. The manuscript reports bulk PBMC RNA-seq, SomaScan proteomics, targeted cytokine measurements, and flow-cytometry profiling, with controlled-access raw data available through EGA under **EGAS50000001393**. This repository is organized so that approved users can retrieve the raw or processed files from EGA and rerun the downstream analyses locally. 

The study design, cohort size, molecular layers, and analytical methods documented here follow the manuscript and supplementary methods: 38 FMF patients were profiled across bulk RNA-seq, cytokine profiling, SomaScan proteomics, and flow cytometry; bulk RNA-seq used kallisto/tximport and DESeq2 with age, sex, and batch/plate covariates when applicable; SomaScan data were normalized with the vendor/Candia workflow and analyzed with limma; targeted cytokines and flow cytometry used nonparametric group comparisons with multiple-testing correction; bulk RNA-seq pathway and TF activities were summarized with clusterProfiler and decoupleR. fileciteturn1file5 fileciteturn1file10

## Repository layout

```text
FMFSub_GitHub_Code/
├── README.md
├── LICENSE
├── environment.yml
├── config/
│   ├── config_template.yaml
│   └── tx2gene_template.tsv
├── metadata/
│   ├── sample_metadata_template.tsv
│   ├── somascan_manifest_template.tsv
│   ├── cytokine_metadata_template.tsv
│   └── flow_metadata_template.tsv
├── scripts/
│   ├── 00_bulkRNA_fastq_to_kallisto.sh
│   ├── 01_bulkRNA_deseq2.R
│   ├── 02_bulkRNA_go_enrichment.R
│   ├── 03_bulkRNA_decoupler.R
│   ├── 04_somascan_limma.R
│   ├── 05_cytokine_stats.R
│   ├── 06_flowcytometry_stats.R
│   └── 07_session_info.R
└── docs/
    └── github_upload_checklist.md
```

## Expected inputs

### 1) Bulk RNA-seq
Use either:
- **raw FASTQ files** downloaded from EGA, followed by `scripts/00_bulkRNA_fastq_to_kallisto.sh`, or
- a **gene-level raw count matrix** already derived from kallisto + tximport.

Required metadata columns in `metadata/sample_metadata.tsv`:
- `sample_id`
- `group` with values `biallelic` or `simple_heterozygous`
- `age`
- `sex`
- optional `batch`
- optional `plate`
- optional additional covariates

### 2) SomaScan proteomics
Provide a tabular matrix with one row per analyte and one column per sample, or one row per sample and one column per analyte. The script can auto-detect orientation if a manifest is supplied. Required metadata columns:
- `sample_id`
- `group`
- optional `age`, `sex`, `plate`, `batch`
- optional annotation columns such as `AptName`, `Target`, `UniProt`

### 3) Targeted cytokines / ELISA / Luminex
Provide a long or wide table. Minimum required fields after reshaping:
- `sample_id`
- `group`
- `analyte`
- `value`

### 4) Flow cytometry
Provide a tabular file containing either:
- frequencies per population, and/or
- MFI values per activation marker and population.

Minimum recommended columns:
- `sample_id`
- `group`
- `feature`
- `value`
- optional `feature_type` (`frequency` or `mfi`)

## Minimal run order

```bash
# 0) create environment
conda env create -f environment.yml
conda activate fmfsub

# 1) copy and edit templates
cp config/config_template.yaml config/config.yaml
cp metadata/sample_metadata_template.tsv metadata/sample_metadata.tsv

# 2) optional raw RNA-seq quantification from FASTQ
bash scripts/00_bulkRNA_fastq_to_kallisto.sh config/config.yaml

# 3) downstream analyses
Rscript scripts/01_bulkRNA_deseq2.R config/config.yaml
Rscript scripts/02_bulkRNA_go_enrichment.R config/config.yaml
Rscript scripts/03_bulkRNA_decoupler.R config/config.yaml
Rscript scripts/04_somascan_limma.R config/config.yaml
Rscript scripts/05_cytokine_stats.R config/config.yaml
Rscript scripts/06_flowcytometry_stats.R config/config.yaml
Rscript scripts/07_session_info.R results/session_info.txt
```

## Main outputs

The scripts write to `results/` and generate:
- differential expression tables for bulk RNA-seq
- VST PCA and sample distance plots
- volcano plots and heatmaps
- GO enrichment tables and dot plots
- decoupleR pathway/TF activity tables and plots
- SomaScan differential abundance tables and volcano plots
- cytokine boxplots with Mann–Whitney tests and BH correction
- flow-cytometry summary tables and plots
- a session-info text file for reproducibility

## Statistical conventions encoded in the scripts

These defaults follow the methods used in the manuscript:
- **bulk RNA-seq**: keep genes with at least 10 counts in at least 20% of samples; DESeq2 Wald test; BH-adjusted FDR; optional `apeglm` LFC shrinkage; significance default `padj < 0.05` with a display threshold `|log2FC| > 1`. fileciteturn1file10
- **SomaScan**: log2-transformed normalized RFU values; limma linear modeling with optional covariates and BH correction. fileciteturn1file6
- **targeted cytokines and flow cytometry**: two-sided Mann–Whitney tests for two-group comparisons with Benjamini–Hochberg correction across tested features, matching the figure-level descriptions in the manuscript. fileciteturn1file17
- **pathway analysis**: clusterProfiler GO enrichment; decoupleR with PROGENy and DoRothEA resources as described in the supplementary methods. fileciteturn1file10

## Notes for GitHub release

1. Do **not** upload controlled-access raw data from EGA.
2. Do upload:
   - code
   - configuration templates
   - sample-sheet templates without patient identifiers
   - a small synthetic example if desired
   - exact software versions
3. In the manuscript, the existing GitHub location is cited as `CrossReaction/Immunaid/tree/main/FMFSub`. If you replace or extend that location, update the code-availability statement accordingly. fileciteturn1file16

## Suggested code-availability wording

> Code for the analysis of the FMF genotype-dosage multi-omics datasets is available in this repository. Controlled-access raw data are available from EGA under accession EGAS50000001393, subject to approval by the ImmunAID data access committee.

## Contact

For questions about the code or data structure, designate the lead bioinformatics contact in the repository description or `README`.
