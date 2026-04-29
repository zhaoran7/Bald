#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
res=/data/sph-zhaor/analysis/bald/res/ld
out=/data/sph-zhaor/analysis/bald/res/pick
traits="bald.qt bald12.bt bald13.bt bald14.bt"

mkdir -p $out

first=1
: > $out/pick.tsv

for t in $traits; do
  f=$res/$t/block.tsv
  [[ -f $f ]] || continue
  if [[ $first -eq 1 ]]; then
    awk 'BEGIN{FS=OFS="	"} NR==1 || $7>1' "$f" > $out/pick.tsv
    first=0
  else
    awk 'BEGIN{FS=OFS="	"} NR>1 && $7>1' "$f" >> $out/pick.tsv
  fi
done
