#!/usr/bin/env bash
set -euo pipefail

tr="bald.qt bald12.bt bald13.bt bald14.bt"
src=/mnt/d/bald/res
ref=/mnt/d/bald/data/raw/37/1kg
out=/mnt/d/bald/res2/ld
plink=plink
r2=0.98
kb=1000
win=999999

mkdir -p "$out/ref"

for c in $(seq 1 22); do
  vcf=$ref/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  idx=$out/ref/chr${c}.rsid.tsv
  [[ -f $idx ]] || bcftools query -f '%ID\t%POS\n' "$vcf" | awk '$1!="."' > "$idx"
done

for t in $tr; do
  in=$src/$t/lead.tsv
  od=$out/$t
  mkdir -p "$od/ld"
  [[ -f $in ]] || continue

  printf "lead_chr\tlead_snp\tlead_bp\tin_1kg\tid1k\tpos1k\n" > "$od/lead.tsv"

  awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3}' "$in" | \
  while IFS=$'\t' read -r chr snp bp; do
    idx=$out/ref/chr${chr}.rsid.tsv
    pos=$(awk -v s="$snp" 'BEGIN{FS="\t"} $1==s{print $2; exit}' "$idx")
    if [[ -n "${pos:-}" ]]; then
      printf "%s\t%s\t%s\t1\t%s\t%s\n" "$chr" "$snp" "$bp" "$snp" "$pos" >> "$od/lead.tsv"
    else
      printf "%s\t%s\t%s\t0\t.\t.\n" "$chr" "$snp" "$bp" >> "$od/lead.tsv"
    fi
  done

  printf "trait\tlead_chr\tlead_snp\tlead_bp\tchr\tpos\tsnp\tR2\n" > "$od/ld.tsv"

  awk 'BEGIN{FS=OFS="\t"} NR>1 && $4==1 {print $1,$2,$3,$5}' "$od/lead.tsv" | \
  while IFS=$'\t' read -r chr snp bp id1k; do
    vcf=$ref/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    tmp=$od/ld/chr${chr}.${snp}

    $plink \
      --vcf "$vcf" \
      --double-id \
      --allow-extra-chr \
      --r2 \
      --ld-snp "$id1k" \
      --ld-window-kb $kb \
      --ld-window $win \
      --ld-window-r2 $r2 \
      --out "$tmp" >/dev/null 2>&1 || continue

    [[ -f ${tmp}.ld ]] || continue

    awk -v t="$t" -v c="$chr" -v s="$snp" -v b="$bp" '
    BEGIN{FS="[ \t]+"; OFS="\t"}
    NR>1 {print t,c,s,b,$5,$6,$7,$8}
    ' "${tmp}.ld" >> "$od/ld.tsv"
  done

  awk '!seen[$0]++' "$od/ld.tsv" > "$od/.tmp" && mv "$od/.tmp" "$od/ld.tsv"

  awk '
    BEGIN{FS=OFS="\t"}
    NR==1{next}
    {
      k=$1 FS $2 FS $3 FS $4
      if(!(k in mn) || $6<mn[k]) mn[k]=$6
      if(!(k in mx) || $6>mx[k]) mx[k]=$6
      u[k FS $7]=1
    }
    END{
      print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"
      for(k in mn){
        n=0
        for(i in u) if(index(i,k FS)==1) n++
        split(k,x,FS)
        print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1
      }
    }' "$od/ld.tsv" | sort -k2,2n -k4,4n > "$od/block.tsv"
done
