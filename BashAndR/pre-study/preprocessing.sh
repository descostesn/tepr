#! /usr/bin/sh

cd /g/romebioinfo/Projects/tepr/downloads/annotations

## Create bed of protein coding genes from gencode file
grep -w transcript gencode.v43.basic.annotation.gtf | grep -w MANE_Select > temp1.gtf

awk 'OFS="\t" {print $0}' MANE_Select.protein_coding.gtf | tr -d '";' | \
sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | \
awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > MANE_Select.protein_coding.bed

rm tmp1.gtf 

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
## NOTE: Getting warning messages -
## Interval chr**:**-** is smaller than the number of windows requested. Skipping.
mkdir makewindow

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' MANE_Select.protein_coding.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v > tmp2.bed

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' MANE_Select.protein_coding.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v  | \
bedtools makewindows -n 200 -i srcwinnum -b stdin > tmp3.bed  

rm tmp2.bed tmp3.bed


awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' MANE_Select.protein_coding.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v  | \
bedtools makewindows -n 200 -i srcwinnum -b stdin | \
sort -k1,1 -k2,2n > makewindow/v43.MANE_protein.window200.bed  

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' Ensembl_canonical_TSL123.lncRNA.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v > tmp4.bed

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' Ensembl_canonical_TSL123.lncRNA.bed | \
bedtools intersect -a stdin -b hg38-blacklist.v2.bed -v  | \
bedtools makewindows -n 200 -i srcwinnum -b stdin | \
sort -k1,1 -k2,2n > makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed 

rm tmp4.bed


## Associating scores to bed
## The umap track was downloaded from https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Unique.Mappability.bb
## and converted to bed with https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed
cd /g/romebioinfo/Projects/tepr/downloads/

ext="bg"
umapk50="annotations/k50.Unique.Mappability.bed" # hg38 mappability windows with mapp k50 > 80%
blacklist="annotations/hg38-blacklist.v2.bed"
window=200

## Processing protein coding
windowS="annotations/makewindow/v43.MANE_protein.window200.bed" #genes with the window
ANNOTATION="annotations/MANE_Select.protein_coding.bed"

mkdir bedgraphs/withzeros
mkdir bedgraphs/mapHigh

for file in bedgraphs/*.$ext ;
do
    filename=$(basename "$file" .$ext);
    echo "starting file :"
    echo $filename;
    if echo $filename | egrep -q "reverse|minus" ;  then
        strand="-"
    elif echo $filename | egrep -q "forward|plus" ; then
        strand='+'
    fi

    echo "removing blacklist region"
    bedtools intersect -a bedgraphs/${filename}.$ext -b <( awk -F "\t" -v OFS="\t" -v myvar=$strand '{if ($6==myvar) print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | \
    bedtools intersect -a stdin -b $blacklist -v) | sort -k1,1 -k2,2n > bedgraphs/withzeros/${filename}.nonzeros.$ext

    echo "removing low mappability region"
    bedtools intersect -a bedgraphs/withzeros/${filename}.nonzeros.${ext} -b $umapk50 -sorted | \
    awk -F "\t" -v OFS="\t" '{print $1,$2,$3,".",$4}' > bedgraphs/mapHigh/${filename}.0.8.$ext

    echo "scoring windows"
    bedmap --echo --wmean --delim "\t" $windowS bedgraphs/mapHigh/${filename}.0.8.$ext | \
    awk -F "_" -v OFS="\t" '{print $1,$2,$3,$4}' | awk -F "\t" -v OFS="\t" -v name="$filename" '{ print $0,$4"_"$5"_"$6"_"$7,name}' | \
    awk -F "\t" -v OFS="\t" '{ print "protein-coding",$1,$2,$3,$4,$5,$6,$7,$9,$10,$8 }' > bedgraphs/${filename}.window${window}.MANE.wmean.name.score ;

    echo "done"
done
rm bedgraphs/withzeros/*
rm bedgraphs/mapHigh/*
mkdir bedgraphs/protein_coding_score
mv bedgraphs/*.score bedgraphs/protein_coding_score/


## Processing lncRNA
windowS="annotations/makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed"
ANNOTATION="annotations/Ensembl_canonical_TSL123.lncRNA.bed"


for file in bedgraphs/*.$ext ;
do
    filename=$(basename "$file" .$ext) ;
    echo "starting file :"
    echo $filename ;

    if echo $filename | egrep -q "reverse|minus" ;  then 
    strand="-"
    elif echo $filename | egrep -q "forward|plus" ; then
    strand='+'
    fi

    echo "removing blacklist region"
    bedtools intersect -a bedgraphs/${filename}.$ext -b <( awk -F "\t" -v OFS="\t" -v myvar=$strand '{if ($6==myvar) print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | \
    bedtools intersect -a stdin -b $blacklist -v) | sort -k1,1 -k2,2n > bedgraphs/withzeros/${filename}.nonzeros.$ext

    echo "removing low mappability region"
    bedtools intersect -a bedgraphs/withzeros/${filename}.nonzeros.$ext -b $umapk50 -sorted | \
    awk -F "\t" -v OFS="\t" '{print $1,$2,$3,".",$4}' > bedgraphs/mapHigh/${filename}.0.8.$ext

    echo "scoring windows"
    bedmap --echo --wmean --delim "\t" $windowS bedgraphs/mapHigh/${filename}.0.8.$ext | \
    awk -F "_" -v OFS="\t" '{print $1,$2,$3,$4}' | awk -F "\t" -v OFS="\t" -v name="$filename" '{ print $0,$4"_"$5"_"$6"_"$7,name}' | \
    awk -F "\t" -v OFS="\t" '{ print "lncRNA",$1,$2,$3,$4,$5,$6,$7,$9,$10,$8 }' > bedgraphs/${filename}.window${window}.MANE.wmean.name.score ;

    echo "done"
done
rm -r bedgraphs/withzeros
rm -r bedgraphs/mapHigh
mkdir bedgraphs/lncRNA_score
mv bedgraphs/*.score bedgraphs/lncRNA_score/





