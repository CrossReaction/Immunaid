#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash scripts/00_bulkRNA_fastq_to_kallisto.sh config/config.yaml"
  exit 1
fi

CONFIG="$1"

get_yaml() {
python - <<PY "$CONFIG" "$1"
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
key = sys.argv[2].split('.')
val = cfg
for k in key:
    val = val[k]
print(val)
PY
}

INDEX=$(get_yaml bulk_rnaseq.fasta_index)
QUANT_DIR=$(get_yaml bulk_rnaseq.quant_dir)
META=$(get_yaml bulk_rnaseq.metadata_file)
SAMPLE_COL=$(get_yaml bulk_rnaseq.sample_id_col)

mkdir -p "$QUANT_DIR"

if [[ ! -f "$INDEX" ]]; then
  echo "kallisto index not found: $INDEX"
  exit 1
fi

if [[ ! -f "$META" ]]; then
  echo "metadata file not found: $META"
  exit 1
fi

echo "Expect metadata to contain: sample_id, fastq_r1, fastq_r2"
echo "Running kallisto quant for each sample listed in $META"

awk 'BEGIN{FS=OFS="\t"} NR==1 {for(i=1;i<=NF;i++){h[$i]=i}; next} {print $h["sample_id"], $h["fastq_r1"], $h["fastq_r2"]}' "$META" | \
while read -r SAMPLE R1 R2; do
  [[ -z "$SAMPLE" ]] && continue
  OUT="$QUANT_DIR/$SAMPLE"
  mkdir -p "$OUT"
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "Skipping $SAMPLE because FASTQ files are missing"
    continue
  fi
  echo "Quantifying $SAMPLE"
  kallisto quant -i "$INDEX" -o "$OUT" -b 100 "$R1" "$R2"
done
