#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
vcf=/data/sph-zhaor/analysis/bald/res/core_vcf
script=$proj/scripts/s4_make_mat.sh
job=$proj/jobs/mat
traits="bald.qt bald12.bt bald13.bt bald14.bt"

mkdir -p $job/cmd $job/out $job/err

for trait in $traits; do
  for chr in $(seq 1 22); do
    n=$(find $vcf/$trait -maxdepth 1 -type d -name "${chr}.*" 2>/dev/null | wc -l)
    [[ $n -gt 0 ]] || continue

    cmd=$job/cmd/${trait}.chr${chr}.sh
    jn=MAT_${trait//./_}_CHR${chr}

    cat > "$cmd" <<EOF
#!/usr/bin/env bash
set -euo pipefail
bash $script $trait $chr
EOF
    chmod +x "$cmd"

    bsub -q ser -n 1 -J $jn       -o $job/out/${trait}.chr${chr}.out       -e $job/err/${trait}.chr${chr}.err       < "$cmd"
  done
done
