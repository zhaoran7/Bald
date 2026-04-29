#!/usr/bin/env bash
set -euo pipefail

root=/data/sph-zhaor/analysis/bald/res/ld
traits="bald.qt bald12.bt bald13.bt bald14.bt"

for trait in $traits; do
  d=$root/$trait
  [[ -d $d ]] || continue

  first=1
  : > "$d/ld.tsv"
  for chr in $(seq 1 22); do
    f=$d/chr${chr}/ld.tsv
    [[ -f $f ]] || continue
    if [[ $first -eq 1 ]]; then
      cat "$f" > "$d/ld.tsv"
      first=0
    else
      awk 'NR>1' "$f" >> "$d/ld.tsv"
    fi
  done
  awk 'NR==1 || $1!="trait"' "$d/ld.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.ld.tmp" && mv "$d/.ld.tmp" "$d/ld.tsv"

  first=1
  : > "$d/block.tsv"
  for chr in $(seq 1 22); do
    f=$d/chr${chr}/block.tsv
    [[ -f $f ]] || continue
    if [[ $first -eq 1 ]]; then
      cat "$f" > "$d/block.tsv"
      first=0
    else
      awk 'NR>1' "$f" >> "$d/block.tsv"
    fi
  done
  awk 'NR==1 || $1!="trait"' "$d/block.tsv" | awk 'NR==1 || !seen[$0]++' > "$d/.block.tmp" && mv "$d/.block.tmp" "$d/block.tsv"
done
