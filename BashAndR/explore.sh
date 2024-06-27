#! /usr/bin/sh

cd /g/romebioinfo/Projects/tepr/downloads/annotations

## Create bed of protein coding genes from gencode file
grep -w transcript gencode.v43.basic.annotation.gtf | grep -w MANE_Select > MANE_Select.protein_coding.gtf

awk 'OFS="\t" {print $0}' MANE_Select.protein_coding.gtf | tr -d '";' | \
sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | \
awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > MANE_Select.protein_coding.bed

rm MANE_Select.protein_coding.gtf

## Create bed of lncRNA genes from gencode file
grep -w transcript gencode.v43.basic.annotation.gtf | grep -w lncRNA | \
grep -w Ensembl_canonical | grep -v not_best_in_genome_evidence | \
grep -v 'transcript_support_level "5"' | \
grep -v 'transcript_support_level "4"' > Ensembl_canonical_TSL123.lncRNA.gtf

awk 'OFS="\t" {print $0}' Ensembl_canonical_TSL123.lncRNA.gtf | tr -d '";' | \
sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | \
awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > Ensembl_canonical_TSL123.lncRNA.bed

rm Ensembl_canonical_TSL123.lncRNA.gtf

## Removing blacklist from bed files. The black list was downloaded from
## https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
mkdir makewindow

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' MANE_Select.protein_coding.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v  | \
bedtools makewindows -n 200 -i srcwinnum -b stdin | \
sort -k1,1 -k2,2n > makewindow/v43.MANE_protein.window200.bed  

