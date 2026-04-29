#!/usr/bin/env bash
set -euo pipefail

trait=$1
chr=$2

proj=/work/sph-zhaor/analysis/bald
res=/data/sph-zhaor/analysis/bald/res
lead=$proj/data/lead/${trait}.lead.tsv
ref=$proj/data/1kg/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
tmp=$res/ld/$trait/chr${chr}
plink=plink

mkdir -p "$tmp"

[[ -f $lead ]] || { echo "missing $lead"; exit 1; }
[[ -f $ref  ]] || { echo "missing $ref";  exit 1; }
which bcftools >/dev/null 2>&1 || { echo "bcftools not found"; exit 1; }

awk -v chr="$chr" 'BEGIN{FS=OFS="	"} NR==1 || $1==chr' "$lead" > "$tmp/lead.tsv"

n=$(awk 'END{print NR-1}' "$tmp/lead.tsv")
[[ $n -gt 0 ]] || {
  printf "trait	lead_chr	lead_snp	lead_bp	chr	pos	snp	R2
" > "$tmp/ld.tsv"
  printf "trait	lead_chr	lead_snp	lead_bp	start	end	n	size_bp
" > "$tmp/block.tsv"
  exit 0
}

printf "trait	lead_chr	lead_snp	lead_bp	chr	pos	snp	R2
" > "$tmp/ld.tsv"

awk 'BEGIN{FS=OFS="	"} NR>1{print $2,$3}' "$tmp/lead.tsv" | while IFS=$'	' read -r snp bp_in; do
  out=$tmp/${snp}

  bp=$(bcftools query -f '%POS
' -i 'ID=="'"$snp"'"' "$ref" | head -n1 || true)
  [[ -n $bp ]] || continue

  printf "%s	%s	%s	%s	%s	%s	%s	1
" "$trait" "$chr" "$snp" "$bp" "$chr" "$bp" "$snp" >> "$tmp/ld.tsv"

  plink     --vcf "$ref"     --double-id     --allow-extra-chr     --r2     --ld-snp "$snp"     --ld-window-kb 1000     --ld-window 999999     --ld-window-r2 0.98     --out "$out" >/dev/null 2>&1 || continue

  [[ -f ${out}.ld ]] || continue

  awk -v t="$trait" -v c="$chr" -v s="$snp" -v b="$bp" '
    BEGIN{FS="[ 	]+"; OFS="	"}
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="CHR_B") ic=i
        else if($i=="BP_B") ip=i
        else if($i=="SNP_B") is=i
        else if($i=="R2") ir=i
      }
      next
    }
    ic && ip && is && ir {print t,c,s,b,$ic,$ip,$is,$ir}
  ' "${out}.ld" >> "$tmp/ld.tsv"
done

awk 'NR==1 || !seen[$0]++' "$tmp/ld.tsv" > "$tmp/.ld.tmp" && mv "$tmp/.ld.tmp" "$tmp/ld.tsv"

awk '
BEGIN{FS=OFS="	"}
NR==1{next}
{
  k=$1 FS $2 FS $3 FS $4
  bp=$4+0
  pos=$6+0
  if(!(k in mn) || bp<mn[k]) mn[k]=bp
  if(!(k in mx) || bp>mx[k]) mx[k]=bp
  if(pos<mn[k]) mn[k]=pos
  if(pos>mx[k]) mx[k]=pos
  u[k FS bp]=1
  u[k FS pos]=1
}
END{
  print "trait","lead_chr","lead_snp","lead_bp","start","end","n","size_bp"
  for(k in mn){
    n=0
    for(i in u) if(index(i,k FS)==1) n++
    split(k,x,FS)
    print x[1],x[2],x[3],x[4],mn[k],mx[k],n,mx[k]-mn[k]+1
  }
}' "$tmp/ld.tsv" | sort -k2,2n -k4,4n > "$tmp/block.tsv"
