#!/usr/bin/env bash
set -euo pipefail

proj=/work/sph-zhaor/analysis/bald
res=/data/sph-zhaor/analysis/bald/res
phy=$res/phy
job=$proj/jobs/phyml
scr=$proj/scripts

mkdir -p $job/cmd $job/out $job/err

which bsub >/dev/null 2>&1 || { echo "bsub not found"; exit 1; }

find $phy -name "*.phy" | sort > $job/phy.tasks.tsv
n=$(wc -l < $job/phy.tasks.tsv)
[[ $n -gt 0 ]] || { echo "no .phy files"; exit 1; }

cmd=$job/cmd/phy.by_region.sh
cat > $cmd <<EOF
#!/usr/bin/env bash
set -euo pipefail
f=\$(sed -n "\${LSB_JOBINDEX}p" $job/phy.tasks.tsv)
phyml -i "\$f" -m HKY85 -c 4 -a e -v e -b 100
EOF
chmod +x $cmd

bsub -q ser -n 1 -J "PHYML[1-$n]" \
  -o $job/out/phyml.%I.out \
  -e $job/err/phyml.%I.err \
  < $cmd
