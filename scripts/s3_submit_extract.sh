#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
pick=/data/sph-zhaor/analysis/bald/res/pick/pick.tsv
script=$proj/scripts/s3_extract_core.sh
job=$proj/jobs/core_vcf
traits="bald.qt bald12.bt bald13.bt bald14.bt"

mkdir -p $job/cmd $job/out $job/err

for trait in $traits; do
  for chr in $(seq 1 22); do
    n=$(awk -v t="$trait" -v c="$chr" 'BEGIN{FS=OFS="	"} NR>1 && $1==t && $2==c {n++} END{print n+0}' "$pick")
    [[ $n -gt 0 ]] || continue

    cmd=$job/cmd/${trait}.chr${chr}.sh
    jn=CORE_${trait//./_}_CHR${chr}

    cat > "$cmd" <<EOF
#!/usr/bin/env bash
set -euo pipefail
bash $script $trait $chr
EOF
    chmod +x "$cmd"

    bsub -q ser -n 1 -J $jn       -o $job/out/${trait}.chr${chr}.out       -e $job/err/${trait}.chr${chr}.err       < "$cmd"
  done
done
