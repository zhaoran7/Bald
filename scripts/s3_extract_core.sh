#!/usr/bin/env bash
set -euo pipefail

trait=$1
chr=$2

proj=/work/sph-zhaor/analysis/bald
pick=/data/sph-zhaor/analysis/bald/res/pick/pick.tsv
out=/data/sph-zhaor/analysis/bald/res/core_vcf
kg=$proj/data/1kg
vindija=$proj/data/Vindija
altai=$proj/data/Altai
chagyr=$proj/data/Chagyr
denisova=$proj/data/Denisova

mkdir -p $out/$trait
which bcftools >/dev/null 2>&1 || { echo "bcftools not found"; exit 1; }

awk -v t="$trait" -v c="$chr" 'BEGIN{FS=OFS="	"} NR>1 && $1==t && $2==c {print $0}' "$pick" | while IFS=$'	' read -r tr lead_chr snp bp st en n size_bp; do
  d=$out/$trait/${lead_chr}.${snp}.${bp}
  lo=$st
  hi=$en
  [[ $bp -lt $lo ]] && lo=$bp
  [[ $bp -gt $hi ]] && hi=$bp
  r=${lead_chr}:${lo}-${hi}
  mkdir -p $d

  bcftools view -r "$r" -m2 -M2 -v snps -Oz -o $d/kg.vcf.gz $kg/ALL.chr${lead_chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  tabix -f -p vcf $d/kg.vcf.gz

  bcftools view -r "$r" -Oz -o $d/vindija.vcf.gz  $vindija/chr${lead_chr}_mq25_mapab100.vcf.gz
  tabix -f -p vcf $d/vindija.vcf.gz

  bcftools view -r "$r" -Oz -o $d/altai.vcf.gz    $altai/chr${lead_chr}_mq25_mapab100.vcf.gz
  tabix -f -p vcf $d/altai.vcf.gz

  bcftools view -r "$r" -Oz -o $d/chagyr.vcf.gz   $chagyr/chr${lead_chr}.noRB.bgz.vcf.gz
  tabix -f -p vcf $d/chagyr.vcf.gz

  bcftools view -r "$r" -Oz -o $d/denisova.vcf.gz $denisova/chr${lead_chr}_mq25_mapab100.vcf.gz
  tabix -f -p vcf $d/denisova.vcf.gz
done
