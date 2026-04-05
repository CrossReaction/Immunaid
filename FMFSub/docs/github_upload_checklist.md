# GitHub upload checklist

## Before pushing
- Replace all templates with finalized non-identifiable metadata templates.
- Verify that no controlled-access EGA raw data or patient identifiers are present.
- Check that `config/config_template.yaml` uses placeholder paths only.
- Run all scripts once on a local dry run.
- Save `sessionInfo()` output and package versions.
- Add repository topics: `FMF`, `multiomics`, `bulk-rnaseq`, `somascan`, `flow-cytometry`, `luminex`.

## Recommended release notes
- cohort: 27 biallelic FMF vs 11 simple heterozygous FMF during acute crisis
- accession: EGAS50000001393
- methods: DESeq2, limma, clusterProfiler, decoupleR, Mann–Whitney + BH correction
- outputs: differential expression, pathway activity, cytokine statistics, flow summary
