#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
script=$proj/scripts/s1_run_ld.sh
job=$proj/jobs/ld
traits="bald.qt bald12.bt bald13.bt bald14.bt"
res=/data/sph-zhaor/analysis/bald/res

mkdir -p $job/cmd $job/out $job/err
mkdir -p $res/ld

for trait in $traits; do
  mkdir -p $job/cmd/$trait $job/out/$trait $job/err/$trait
  mkdir -p $res/ld/$trait

  for chr in {1..22}; do
    cmd=$job/cmd/$trait/chr${chr}.sh
    jn=LD_${trait//./_}_CHR${chr}

    cat > "$cmd" <<EOF
#!/usr/bin/env bash
set -euo pipefail
bash $script $trait $chr
EOF
    chmod +x "$cmd"

    bsub -q ser -n 1 -J $jn       -o $job/out/$trait/chr${chr}.out       -e $job/err/$trait/chr${chr}.err       < "$cmd"
  done
done
