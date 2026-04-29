#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
res=/data/sph-zhaor/analysis/bald/res
mat=$res/mat
scr=$proj/scripts
job=$proj/jobs/phyml

mkdir -p $job/cmd $job/out $job/err
which bsub >/dev/null 2>&1 || { echo "bsub not found"; exit 1; }
which Rscript >/dev/null 2>&1 || { echo "Rscript not found"; exit 1; }

find $mat -mindepth 1 -maxdepth 1 -type d | sort | awk -F/ '{print $NF}' > $job/s5.tasks.tsv
n=$(wc -l < $job/s5.tasks.tsv)
[[ $n -gt 0 ]] || { echo "no s5 trait tasks"; exit 1; }

cmd=$job/cmd/s5.by_trait.sh
cat > $cmd <<EOF
#!/usr/bin/env bash
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r
trait=\$(sed -n "\${LSB_JOBINDEX}p" $job/s5.tasks.tsv)
Rscript $scr/s5_make_hap.R "\$trait"
EOF
chmod +x $cmd

jid=$(bsub -q ser -n 1 -J "S5HAP[1-$n]" \
  -o $job/out/s5.%I.out \
  -e $job/err/s5.%I.err \
  < $cmd | awk -F'[<>]' 'NF>=3{print $2}')
[[ -n $jid ]] || { echo "failed to submit s5 array"; exit 1; }

cmd=$job/cmd/s5.merge.s6.sh
cat > $cmd <<EOF
#!/usr/bin/env bash
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r
Rscript $scr/s5_make_hap.R merge
Rscript $scr/s6_make_phy.R
EOF
chmod +x $cmd

bsub -q ser -n 1 -w "done($jid)" -J S5MERGE \
  -o $job/out/s5_merge.out \
  -e $job/err/s5_merge.err \
  < $cmd
