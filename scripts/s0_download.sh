#!/usr/bin/env bash
set -euo pipefail

root=/mnt/d/bald/data/raw/37
onekg=$root/1kg
vindija=$root/Vindija
altai=$root/Altai
chagyr=$root/Chagyr
denisova=$root/Denisova

mkdir -p "$onekg" "$vindija" "$altai" "$chagyr" "$denisova"

get_one () {
  local url=$1
  local out=$2
  local tmp=${out}.part
  rm -f "$tmp"
  wget -O "$tmp" "$url"
  mv "$tmp" "$out"
}

ok_gz () {
  local f=$1
  [[ -s "$f" ]] || return 1
  gzip -t "$f" >/dev/null 2>&1 || return 1
  return 0
}

ok_tbi () {
  local f=$1
  [[ -s "$f" ]] || return 1
  return 0
}

ok_vcf_bgz () {
  local vcf=$1
  local tbi=$2
  [[ -s "$vcf" ]] || return 1
  [[ -s "$tbi" ]] || return 1
  gzip -t "$vcf" >/dev/null 2>&1 || return 1
  tabix -l "$vcf" >/dev/null 2>&1 || return 1
  return 0
}

echo "[1/6] 1000G v5a"
for c in $(seq 1 22); do
  vcf="$onekg/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  tbi="${vcf}.tbi"
  if ok_vcf_bgz "$vcf" "$tbi"; then
    echo "[OK] 1kg chr${c}"
  else
    echo "[GET] 1kg chr${c}"
    rm -f "$vcf" "$tbi"
    get_one \
      "https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
      "$vcf"
    get_one \
      "https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi" \
      "$tbi"
    ok_vcf_bgz "$vcf" "$tbi"
  fi
done

panel="$onekg/integrated_call_samples_v3.20130502.ALL.panel"
if [[ -s "$panel" ]]; then
  echo "[OK] 1kg panel"
else
  echo "[GET] 1kg panel"
  rm -f "$panel"
  get_one \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" \
    "$panel"
fi

echo "[2/6] Vindija"
for c in $(seq 1 22); do
  vcf="$vindija/chr${c}_mq25_mapab100.vcf.gz"
  tbi="${vcf}.tbi"
  if ok_vcf_bgz "$vcf" "$tbi"; then
    echo "[OK] Vindija chr${c}"
  else
    echo "[GET] Vindija chr${c}"
    rm -f "$vcf" "$tbi"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz" \
      "$vcf"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${c}_mq25_mapab100.vcf.gz.tbi" \
      "$tbi"
    ok_vcf_bgz "$vcf" "$tbi"
  fi
done

echo "[3/6] Altai"
for c in $(seq 1 22); do
  vcf="$altai/chr${c}_mq25_mapab100.vcf.gz"
  tbi="${vcf}.tbi"
  if ok_vcf_bgz "$vcf" "$tbi"; then
    echo "[OK] Altai chr${c}"
  else
    echo "[GET] Altai chr${c}"
    rm -f "$vcf" "$tbi"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz" \
      "$vcf"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${c}_mq25_mapab100.vcf.gz.tbi" \
      "$tbi"
    ok_vcf_bgz "$vcf" "$tbi"
  fi
done

echo "[4/6] Chagyr"
for c in $(seq 1 22); do
  raw="$chagyr/chr${c}.noRB.vcf.gz"
  raw_tbi="${raw}.tbi"
  bgz="$chagyr/chr${c}.noRB.bgz.vcf.gz"
  bgz_tbi="${bgz}.tbi"

  if ok_vcf_bgz "$bgz" "$bgz_tbi"; then
    echo "[OK] Chagyr chr${c}"
  else
    echo "[GET] Chagyr chr${c}"
    rm -f "$raw" "$raw_tbi" "$bgz" "$bgz_tbi"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz" \
      "$raw"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz.tbi" \
      "$raw_tbi"
    ok_gz "$raw"
    gunzip -c "$raw" | bgzip -c > "$bgz"
    tabix -f -p vcf "$bgz"
    ok_vcf_bgz "$bgz" "$bgz_tbi"
  fi
done

echo "[5/6] Denisova"
for c in $(seq 1 22); do
  vcf="$denisova/chr${c}_mq25_mapab100.vcf.gz"
  tbi="${vcf}.tbi"
  if ok_vcf_bgz "$vcf" "$tbi"; then
    echo "[OK] Denisova chr${c}"
  else
    echo "[GET] Denisova chr${c}"
    rm -f "$vcf" "$tbi"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz" \
      "$vcf"
    get_one \
      "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${c}_mq25_mapab100.vcf.gz.tbi" \
      "$tbi"
    ok_vcf_bgz "$vcf" "$tbi"
  fi
done

echo "[6/6] quick check chr3 lead SNP"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
  "$onekg/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
| awk '$3=="rs35044562"{print "[OK] rs35044562\t"$0; found=1} END{if(!found) print "[WARN] rs35044562 not found"}'

echo "[DONE]"
