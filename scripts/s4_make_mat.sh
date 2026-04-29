#!/usr/bin/env bash
set -euo pipefail

trait=$1
chr=$2

vcf=/data/sph-zhaor/analysis/bald/res/core_vcf/$trait
mat=/data/sph-zhaor/analysis/bald/res/mat/$trait

mkdir -p $mat
which bcftools >/dev/null 2>&1 || { echo "bcftools not found"; exit 1; }

for d in $vcf/${chr}.*; do
  [[ -d $d ]] || continue
  id=$(basename $d)
  o=$mat/$id
  mkdir -p $o

  bcftools query -f '%CHROM	%POS	%REF	%ALT	%INFO/AA[	%GT]
' $d/kg.vcf.gz > $o/kg.tsv
  bcftools query -f '%CHROM	%POS	%REF	%ALT[	%GT]
' $d/vindija.vcf.gz > $o/vindija.tsv
  bcftools query -f '%CHROM	%POS	%REF	%ALT[	%GT]
' $d/altai.vcf.gz > $o/altai.tsv
  bcftools query -f '%CHROM	%POS	%REF	%ALT[	%GT]
' $d/chagyr.vcf.gz > $o/chagyr.tsv
  bcftools query -f '%CHROM	%POS	%REF	%ALT[	%GT]
' $d/denisova.vcf.gz > $o/denisova.tsv
done
