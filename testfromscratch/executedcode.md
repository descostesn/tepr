The first piece of code was executed from `testfromscratch/` in which the annotation file `gencode.v43.basic.annotation.gtf` was copied. The folder hence contains only one file at this point. Retrieving the code from `BashAndR\pre-study\documentation\TSV_creation_lite_20240617.html`, and replacing `*gtf` by `gencode.v43.basic.annotation.gtf` (otherwise it generates an infinite loop), the first piece of code executed is:

```
#!/usr/bin/sh

grep -w transcript gencode.v43.basic.annotation.gtf | grep -w MANE_Select > MANE_Select.protein_coding.gtf

awk 'OFS="\t" {print $0}' MANE_Select.protein_coding.gtf | tr -d '";' | sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > MANE_Select.protein_coding.bed
```

Using the same file, retrieving the transcripts:

```
#!/usr/bin/sh

grep -w transcript gencode.v43.basic.annotation.gtf | grep -w lncRNA | grep -w Ensembl_canonical | grep -v not_best_in_genome_evidence | grep -v 'transcript_support_level "5"' | grep -v 'transcript_support_level "4"' > Ensembl_canonical_TSL123.lncRNA.gtf

awk 'OFS="\t" {print $0}' Ensembl_canonical_TSL123.lncRNA.gtf | tr -d '";' | sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > Ensembl_canonical_TSL123.lncRNA.bed
```

The files produced have the following number of lines:

```
> wc -l gencode.v43.basic.annotation.gtf MANE_Select.protein_coding.gtf MANE_Select.protein_coding.bed Ensembl_canonical_TSL123.lncRNA.gtf Ensembl_canonical_TSL123.lncRNA.bed

1998585 gencode.v43.basic.annotation.gtf
19077 MANE_Select.protein_coding.gtf
19077 MANE_Select.protein_coding.bed
14271 Ensembl_canonical_TSL123.lncRNA.gtf
14271 Ensembl_canonical_TSL123.lncRNA.bed
```

In the code below, note that the file of Victor is called `hg38-blacklist.v2.sorted.bed` and that therefore an extra step of sorting might have been performed. The awk command throws a warning `WARNING: Interval chrX:135309480-135309659 is smaller than the number of windows requested. Skipping.`. The echo command returns `./makewindow/v43.MANE_protein.window200.bed`. Removing the black list and make windows from the protein coding file:

```
#!/usr/bin/sh

ANNOTATION="MANE_Select.protein_coding.bed"
WORKING="."
blacklist="hg38-blacklist.v2.bed"
window=200

mkdir makewindow
awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | bedtools intersect -a stdin -b $blacklist -v  | bedtools makewindows -n $window -i srcwinnum -b stdin | sort -k1,1 -k2,2n > $WORKING/makewindow/v43.MANE_protein.window${window}.bed  

echo $WORKING/makewindow/v43.MANE_protein.window${window}.bed 
```

Performing the same thing for lncRNA. The echo command returns `./makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed`:

```
#!/usr/bin/sh

ANNOTATION="Ensembl_canonical_TSL123.lncRNA.bed"
WORKING="."
blacklist="hg38-blacklist.v2.bed"
window=200

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | bedtools intersect -a stdin -b $blacklist -v  | bedtools makewindows -n $window -i srcwinnum -b stdin | sort -k1,1 -k2,2n > $WORKING/makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed 

echo $WORKING/makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed
```

Note that the awk command above gives the following warnings:

```
WARNING: Interval chr1:111181374-111181491 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:16905199-16905396 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:17439186-17439328 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:222742640-222742729 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:247189576-247189725 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:27773858-27774041 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr1:28648600-28648730 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr10:13099855-13099968 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr10:133374736-133374869 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr11:123741776-123741896 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr11:65492897-65493058 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr11:94650324-94650442 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr12:121874193-121874337 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr12:4275562-4275755 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr13:113165002-113165183 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr15:40874433-40874595 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr15:58456188-58456348 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr16:173502-173660 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr16:177313-177490 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr16:29709156-29709328 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr16:73014319-73014476 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:19868568-19868689 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:2748078-2748182 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:30117079-30117172 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:41867581-41867736 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:57955965-57956143 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:58659125-58659306 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr17:80149627-80149798 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr18:8406761-8406953 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr2:25362232-25362427 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr2:61854376-61854552 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr2:64275361-64275503 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr2:85387074-85387146 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr20:38961925-38962111 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr22:38736610-38736792 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr3:179584271-179584410 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr3:194637505-194637664 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr3:58329965-58330118 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr4:158841813-158841980 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr4:68376551-68376665 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr5:139775305-139775472 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:156780327-156780455 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:26474337-26474441 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:28335662-28335828 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:35764401-35764503 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:41764292-41764460 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:41791410-41791477 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:43770429-43770616 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr6:81969453-81969539 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr7:148941483-148941638 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr7:16471184-16471373 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr7:25855322-25855483 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr9:10948372-10948481 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr9:129328261-129328401 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr9:4850299-4850373 is smaller than the number of windows requested. Skipping.
WARNING: Interval chr9:97214834-97214929 is smaller than the number of windows requested. Skipping.
```

The produced files have the follwing number of lines

```
> wc -l makewindow/v43.MANE_protein.window200.bed makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed
3725600 makewindow/v43.MANE_protein.window200.bed
2772000 makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed
```

Now considering the following bedgraph files that were copied to the folder `./bedgraph255`, numbers of lines are indicated below:

```
14857004 ctrl_rep1.forward.bg
14095382 ctrl_rep1.reverse.bg
17346714 ctrl_rep2.forward.bg
16515469 ctrl_rep2.reverse.bg
11037780 HS_rep1.forward.bg
10605325 HS_rep1.reverse.bg
16466636 HS_rep2.forward.bg
15766282 HS_rep2.reverse.bg
```

The following code that was copied in a file `scoring.zsh` and is dealing with protein coding genes:

```
#!/usr/bin/zsh

## DO NOT FORGET TO CHANGE THE EXTENSION OF THE BEDGRAPH IF NEEDED
ext="bg"
##it works and validated by looking at head and tail of documents and exemple genes
WORKING="./bedgraph255"
window=200
######
umapk50="k50.umap.hg38.0.8.bed" # hg38 mappability windows with mapp k50 > 80%
windowS="makewindow/v43.MANE_protein.window${window}.bed" #genes with the window
blacklist="hg38-blacklist.v2.bed"
ANNOTATION="MANE_Select.protein_coding.bed"
#######

mkdir $WORKING/withzeros
mkdir $WORKING/mapHigh

for file in $WORKING/*.$ext ;
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
 bedtools intersect \
-a $WORKING/${filename}.$ext \
-b <( awk -F "\t" -v OFS="\t" -v myvar=$strand '{if ($6==myvar) print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | bedtools intersect -a stdin -b $blacklist -v) \
 | sort -k1,1 -k2,2n > $WORKING/withzeros/${filename}.nonzeros.$ext


echo "removing low mappability region"


bedtools intersect -a $WORKING/withzeros/${filename}.nonzeros.${ext} -b $umapk50 -sorted | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,".",$4}' > $WORKING/mapHigh/${filename}.0.8.$ext

echo "scoring windows"

bedmap --echo --wmean --delim "\t" $windowS $WORKING/mapHigh/${filename}.0.8.$ext | awk -F "_" -v OFS="\t" '{print $1,$2,$3,$4}' | awk -F "\t" -v OFS="\t" -v name="$filename" '{ print $0,$4"_"$5"_"$6"_"$7,name}' | awk -F "\t" -v OFS="\t" '{ print "protein-coding",$1,$2,$3,$4,$5,$6,$7,$9,$10,$8 }' > $WORKING/${filename}.window${window}.MANE.wmean.name.score ;

echo "done"

done

mkdir $WORKING/protein_coding_score
mv $WORKING/*.score $WORKING/protein_coding_score/
mv $WORKING/withzeros $WORKING/withzeros-proteincoding
mv $WORKING/mapHigh $WORKING/mapHigh-proteincoding
```

The code above was executed with the command:

```
#!/usr/bin/sh

zsh scoring.zsh
```

The command gave the output:

```
starting file :
ctrl_rep1.forward
removing blacklist region
***** WARNING: File ./bedgraph255/ctrl_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       80      1.0061

***** WARNING: File ./bedgraph255/ctrl_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       80      1.0061

removing low mappability region
scoring windows
done
starting file :
ctrl_rep1.reverse
removing blacklist region
***** WARNING: File ./bedgraph255/ctrl_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       580     0

***** WARNING: File ./bedgraph255/ctrl_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       580     0

removing low mappability region
scoring windows
done
starting file :
ctrl_rep2.forward
removing blacklist region
***** WARNING: File ./bedgraph255/ctrl_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       470     0

***** WARNING: File ./bedgraph255/ctrl_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       470     0

removing low mappability region
scoring windows
done
starting file :
ctrl_rep2.reverse
removing blacklist region
***** WARNING: File ./bedgraph255/ctrl_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       610     0

***** WARNING: File ./bedgraph255/ctrl_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       610     0

removing low mappability region
scoring windows
done
starting file :
HS_rep1.forward
removing blacklist region
***** WARNING: File ./bedgraph255/HS_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       760     0

***** WARNING: File ./bedgraph255/HS_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       760     0

removing low mappability region
scoring windows
done
starting file :
HS_rep1.reverse
removing blacklist region
***** WARNING: File ./bedgraph255/HS_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       2290    0

***** WARNING: File ./bedgraph255/HS_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       2290    0

removing low mappability region
scoring windows
done
starting file :
HS_rep2.forward
removing blacklist region
***** WARNING: File ./bedgraph255/HS_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       840     0

***** WARNING: File ./bedgraph255/HS_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       840     0

removing low mappability region
scoring windows
done
starting file :
HS_rep2.reverse
removing blacklist region
***** WARNING: File ./bedgraph255/HS_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       40      0

***** WARNING: File ./bedgraph255/HS_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       40      0

removing low mappability region
scoring windows
done
```

The files produced have the following numbers of lines:

```
> wc -l bedgraph255/withzeros-proteincoding/*
  11833015 bedgraph255/withzeros-proteincoding/ctrl_rep1.forward.nonzeros.bg
  11028069 bedgraph255/withzeros-proteincoding/ctrl_rep1.reverse.nonzeros.bg
  13764078 bedgraph255/withzeros-proteincoding/ctrl_rep2.forward.nonzeros.bg
  12867337 bedgraph255/withzeros-proteincoding/ctrl_rep2.reverse.nonzeros.bg
   8299043 bedgraph255/withzeros-proteincoding/HS_rep1.forward.nonzeros.bg
   7801098 bedgraph255/withzeros-proteincoding/HS_rep1.reverse.nonzeros.bg
  12597790 bedgraph255/withzeros-proteincoding/HS_rep2.forward.nonzeros.bg
  11807275 bedgraph255/withzeros-proteincoding/HS_rep2.reverse.nonzeros.bg

> wc -l bedgraph255/mapHigh-proteincoding/*
  20651207 bedgraph255/mapHigh-proteincoding/ctrl_rep1.forward.0.8.bg
  19098231 bedgraph255/mapHigh-proteincoding/ctrl_rep1.reverse.0.8.bg
  22434617 bedgraph255/mapHigh-proteincoding/ctrl_rep2.forward.0.8.bg
  20800178 bedgraph255/mapHigh-proteincoding/ctrl_rep2.reverse.0.8.bg
  17325729 bedgraph255/mapHigh-proteincoding/HS_rep1.forward.0.8.bg
  16056740 bedgraph255/mapHigh-proteincoding/HS_rep1.reverse.0.8.bg
  21345438 bedgraph255/mapHigh-proteincoding/HS_rep2.forward.0.8.bg
  19808651 bedgraph255/mapHigh-proteincoding/HS_rep2.reverse.0.8.bg

> wc -l bedgraph255/protein_coding_score/*
   3725600 bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/ctrl_rep1.reverse.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/ctrl_rep2.forward.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/ctrl_rep2.reverse.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/HS_rep1.forward.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/HS_rep1.reverse.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/HS_rep2.forward.window200.MANE.wmean.name.score
   3725600 bedgraph255/protein_coding_score/HS_rep2.reverse.window200.MANE.wmean.name.score
```

The following code that was copied in a file `lncRNA.zsh` and is dealing with lncRNAs:

```
#!/usr/bin/zsh

ext="bg"
WORKING="bedgraph255"
window=200

umapk50="k50.umap.hg38.0.8.bed"
windowS="makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed"
blacklist="hg38-blacklist.v2.bed"
ANNOTATION="Ensembl_canonical_TSL123.lncRNA.bed"

mkdir $WORKING/withzeros
mkdir $WORKING/mapHigh

for file in $WORKING/*.$ext ;
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
 bedtools intersect \
-a $WORKING/${filename}.$ext \
-b <( awk -F "\t" -v OFS="\t" -v myvar=$strand '{if ($6==myvar) print $1,$2,$3,$4"_"$5"_"$6}' $ANNOTATION | bedtools intersect -a stdin -b $blacklist -v) \
 | sort -k1,1 -k2,2n > $WORKING/withzeros/${filename}.nonzeros.$ext


echo "removing low mappability region"

bedtools intersect -a $WORKING/withzeros/${filename}.nonzeros.$ext -b $umapk50 -sorted | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,".",$4}' > $WORKING/mapHigh/${filename}.0.8.$ext

echo "scoring windows"

bedmap --echo --wmean --delim "\t" $windowS $WORKING/mapHigh/${filename}.0.8.$ext | awk -F "_" -v OFS="\t" '{print $1,$2,$3,$4}' | awk -F "\t" -v OFS="\t" -v name="$filename" '{ print $0,$4"_"$5"_"$6"_"$7,name}' | awk -F "\t" -v OFS="\t" '{ print "lncRNA",$1,$2,$3,$4,$5,$6,$7,$9,$10,$8 }' > $WORKING/${filename}.window${window}.MANE.wmean.name.score ;

echo "done"

done

mkdir $WORKING/lncRNA_score
mv $WORKING/*.score $WORKING/lncRNA_score/
mv $WORKING/withzeros $WORKING/withzeros-lncRNA
mv $WORKING/mapHigh $WORKING/mapHigh-lncRNA
```

The code above was executed with the command:

```
#!/usr/bin/sh

zsh lncRNA.zsh
```

The command gave the output:

```
starting file :
ctrl_rep1.forward
removing blacklist region
***** WARNING: File bedgraph255/ctrl_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       80      1.0061

***** WARNING: File bedgraph255/ctrl_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       80      1.0061

removing low mappability region
scoring windows
done
starting file :
ctrl_rep1.reverse
removing blacklist region
***** WARNING: File bedgraph255/ctrl_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       580     0

***** WARNING: File bedgraph255/ctrl_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       580     0

removing low mappability region
scoring windows
done
starting file :
ctrl_rep2.forward
removing blacklist region
***** WARNING: File bedgraph255/ctrl_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       470     0

***** WARNING: File bedgraph255/ctrl_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       470     0

removing low mappability region
scoring windows
done
starting file :
ctrl_rep2.reverse
removing blacklist region
***** WARNING: File bedgraph255/ctrl_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       610     0

***** WARNING: File bedgraph255/ctrl_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       610     0

removing low mappability region
scoring windows
done
starting file :
HS_rep1.forward
removing blacklist region
***** WARNING: File bedgraph255/HS_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       760     0

***** WARNING: File bedgraph255/HS_rep1.forward.bg has inconsistent naming convention for record:
GL000008.2      0       760     0

removing low mappability region
scoring windows
done
starting file :
HS_rep1.reverse
removing blacklist region
***** WARNING: File bedgraph255/HS_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       2290    0

***** WARNING: File bedgraph255/HS_rep1.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       2290    0

removing low mappability region
scoring windows
done
starting file :
HS_rep2.forward
removing blacklist region
***** WARNING: File bedgraph255/HS_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       840     0

***** WARNING: File bedgraph255/HS_rep2.forward.bg has inconsistent naming convention for record:
GL000008.2      0       840     0

removing low mappability region
scoring windows
done
starting file :
HS_rep2.reverse
removing blacklist region
***** WARNING: File bedgraph255/HS_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       40      0

***** WARNING: File bedgraph255/HS_rep2.reverse.bg has inconsistent naming convention for record:
GL000008.2      0       40      0

removing low mappability region
scoring windows
done
```

The files produced have the following numbers of lines:

```
> wc -l bedgraph255/withzeros-lncRNA/*
   426766 bedgraph255/withzeros-lncRNA/ctrl_rep1.forward.nonzeros.bg
   431670 bedgraph255/withzeros-lncRNA/ctrl_rep1.reverse.nonzeros.bg
   495350 bedgraph255/withzeros-lncRNA/ctrl_rep2.forward.nonzeros.bg
   514695 bedgraph255/withzeros-lncRNA/ctrl_rep2.reverse.nonzeros.bg
   353893 bedgraph255/withzeros-lncRNA/HS_rep1.forward.nonzeros.bg
   358254 bedgraph255/withzeros-lncRNA/HS_rep1.reverse.nonzeros.bg
   483067 bedgraph255/withzeros-lncRNA/HS_rep2.forward.nonzeros.bg
   497714 bedgraph255/withzeros-lncRNA/HS_rep2.reverse.nonzeros.bg

> wc -l bedgraph255/mapHigh-lncRNA/*
  3716631 bedgraph255/mapHigh-lncRNA/ctrl_rep1.forward.0.8.bg
  3421650 bedgraph255/mapHigh-lncRNA/ctrl_rep1.reverse.0.8.bg
  3777258 bedgraph255/mapHigh-lncRNA/ctrl_rep2.forward.0.8.bg
  3494694 bedgraph255/mapHigh-lncRNA/ctrl_rep2.reverse.0.8.bg
  3649638 bedgraph255/mapHigh-lncRNA/HS_rep1.forward.0.8.bg
  3355952 bedgraph255/mapHigh-lncRNA/HS_rep1.reverse.0.8.bg
  3763918 bedgraph255/mapHigh-lncRNA/HS_rep2.forward.0.8.bg
  3477181 bedgraph255/mapHigh-lncRNA/HS_rep2.reverse.0.8.bg

> wc -l bedgraph255/lncRNA_score/*
   2772000 bedgraph255/lncRNA_score/ctrl_rep1.forward.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/ctrl_rep1.reverse.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/ctrl_rep2.forward.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/ctrl_rep2.reverse.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/HS_rep1.forward.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/HS_rep1.reverse.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/HS_rep2.forward.window200.MANE.wmean.name.score
   2772000 bedgraph255/lncRNA_score/HS_rep2.reverse.window200.MANE.wmean.name.score
```

Executing the following R code to combine the files for protein coding:

```
main_directory <- "bedgraph255"
write_file_protein_coding <- paste0(main_directory,"/protCoding_dTAG_Cugusi_stranded_20230810.tsv")
write_file_lncRNA <- paste0(main_directory,"/lncRNA_dTAG_Cugusi_stranded_20230810.tsv")
Big_tsv <-paste0(main_directory,"/dTAG_Cugusi_stranded_20230810.tsv") 

library(dplyr)
library(purrr)

### DO NOT FORGET TO CHANGE THE BEDGRAPH EXTENSION TO TH ONE IN USE IN THE FOLDER
window <- 200
working_directory <- paste0(main_directory,"/protein_coding_score")
bedgraph_files <- list.files(main_directory, pattern = "*.bg", full.names = TRUE)

files <- bedgraph_files %>%
  map(~{
    filename <- tools::file_path_sans_ext(basename(.))
    file.path(working_directory, paste0(filename, ".window", window, ".MANE.wmean.name.score"))
  })

# reading all the files
list_of_dfs <- lapply(files, read.delim, header= F, sep= "\t",  na.strings = "NAN", dec = ".", col.names = c("biotype","chr", "coor1", "coor2","transcript", "gene", "strand","window","id","dataset","score"), stringsAsFactors=FALSE)

# joining all the files
joined_df <- purrr::reduce(list_of_dfs, dplyr::left_join, by = c("biotype","chr","coor1","coor2","transcript","gene","strand","window","id")) %>% 
   dplyr::filter(strand!="Y") ## the last filter remove the PAR genes (pseudoautosomal genes both in X and Y)

 write.table(joined_df, file = write_file_protein_coding, sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)

 rm(list_of_dfs)
 rm(joined_df)
```

The code above was copied to the file `joinprotcod.R` and executed with R-4.4.1:

```
Rscript joinprotcod.R
```

The following file was obtained:

```
> wc -l bedgraph255/protCoding_dTAG_Cugusi_stranded_20230810.tsv
3722600 bedgraph255/protCoding_dTAG_Cugusi_stranded_20230810.tsv
```

Executing the following R code to combine the files for lncRNA:

```
library(dplyr)
library(purrr)

### DO NOT FORGET TO CHANGE THE BEDGRAPH EXTENSION TO TH ONE IN USE IN THE FOLDER
window <- 200
main_directory <- "bedgraph255"
working_directory <- paste0(main_directory,"/lncRNA_score")
write_file_lncRNA <- paste0(main_directory,"/lncRNA_dTAG_Cugusi_stranded_20230810.tsv")

bedgraph_files <- list.files(main_directory, pattern = "*.bg", full.names = TRUE)

files <- bedgraph_files %>%
  map(~{
    filename <- tools::file_path_sans_ext(basename(.))
    file.path(working_directory, paste0(filename, ".window", window, ".MANE.wmean.name.score"))
  })
# reading all the files
list_of_dfs <- lapply(files, read.delim, header= F, sep= "\t",  na.strings = "NAN", dec = ".", col.names = c("biotype","chr", "coor1", "coor2","transcript", "gene", "strand","window","id","dataset","score"), stringsAsFactors=FALSE)
# joining all the files
joined_df <- purrr::reduce(list_of_dfs, dplyr::left_join, by = c("biotype","chr","coor1","coor2","transcript","gene","strand","window","id")) %>% 
   dplyr::filter(strand!="Y") ## the last filter remove the PAR genes (pseudoautosomal genes both in X and Y)

write.table(joined_df, file = write_file_lncRNA, sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
rm(list_of_dfs)
rm(joined_df)
```

The code above was copied to the file `joinlncrna.R` and executed with R-4.4.1:

```
Rscript joinlncrna.R
```

The following file was obtained:

```
> wc -l bedgraph255/lncRNA_dTAG_Cugusi_stranded_20230810.tsv
    2770600 bedgraph255/lncRNA_dTAG_Cugusi_stranded_20230810.tsv
```

The R code below was used to merge the tables:

```
library(dplyr)
library(purrr)

main_directory <- "bedgraph255"
Big_tsv <-paste0(main_directory,"/dTAG_Cugusi_stranded_20230810.tsv") 

tsv_files <- list.files(main_directory, pattern = "*.tsv", full.names = TRUE)

files <- tsv_files %>%
  map(~{
    filename <- tools::file_path_sans_ext(basename(.))
    file.path(main_directory, paste0(filename, ".tsv"))
  })
# reading all the files
list_of_dfs <- lapply(files, read.delim, header= F, sep= "\t",  na.strings = "NAN", dec = ".", stringsAsFactors=FALSE)
# joining all the files
bound_df <- purrr::reduce(list_of_dfs, dplyr::bind_rows)

write.table(bound_df, file = Big_tsv, sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
rm(list_of_dfs)
rm(bound_df)
```

The code was copied to `jointsvs.R` and executed with the command:

```
Rscript jointsvs.R
```

The following file was obtained:

```
> wc -l bedgraph255/dTAG_Cugusi_stranded_20230810.tsv
    6493200 bedgraph255/dTAG_Cugusi_stranded_20230810.tsv
```


Note that in the provided script, the input table has a different name: `Cugusi_protein-lncRNA_stranded_analysis_MAPQ255_20230705.tsv`. It was replaced by the one obtained above `dTAG_Cugusi_stranded_20230810.tsv` in the variable `name_table`. The following code reads and filters this file based on a minimal expression threshold:

```
#setting columns names
#files names have to be : condition_rep#.strand, a condition can contain several sub condition (cond1a_cond1b)

working_directory <- "bedgraph255" 
## in the working directory the files are: 
# ctrl_rep1.forward.bg  ctrl_rep2.forward.bg  HS_rep1.forward.bg  HS_rep2.forward.bg
# ctrl_rep1.reverse.bg  ctrl_rep2.reverse.bg  HS_rep1.reverse.bg  HS_rep2.reverse.bg
# Cugusi_protein-lncRNA_stranded_analysis_MAPQ255_20230705.tsv 

extension <- "*.bg"
name_table <- "dTAG_Cugusi_stranded_20230810.tsv"
rounding <- 10

getting_var_names <- function(extension, working_directory) {
  
# This function uses the extension and working directory to get the condition names, the number of replicates, and the variable names.
# It needs the file names to be written in the form:
# Condition_rep#.strand.extension such as :
# HS_rep1.reverse.bg
  
  
# In input: extension such as "*.bg" and the working directory
  
  
  bedgraph_files <- list.files(working_directory, pattern = extension, full.names = TRUE)
  files <- bedgraph_files %>%
    map(~{
      filename <- tools::file_path_sans_ext(basename(.))
    })
  
  string <- files
  var_names <- string
  
  for (i in seq_along(var_names)) {
    if (grepl("(reverse|forward)", var_names[i])) {
      var_names[i] <- gsub("reverse", "minus", var_names[i])
      var_names[i] <- gsub("forward", "plus", var_names[i])
    }
  }
  
  # Extract conditions
  Conditions <- unique(sub("(\\w+)_rep\\d+.*", "\\1", var_names)) ## verify it can work with several "_" 
  
  # Extract replication numbers
  replicate_numbers <- unique(sub(".*_rep(\\d+).*", "\\1", var_names))
  
  ###########
  # Setting fixed column names
  fixed_names <- c("biotype","chr", "coor1", "coor2","transcript", "gene", "strand","window","id") #biotype is added in the .tsv file
  
  num_fixed_cols <- length(fixed_names)
  # Number of variable columns based on input file
  num_var_cols <- ncol(table) - num_fixed_cols
  #Generate variable column names for category 2, alternating with category 1
  test <- rep(var_names, each=2)
  suffix <- rep(c("", "_score"), length(var_names))
  test <- paste0(test, suffix)
  col_names <- c(fixed_names, test)
  
  return(list(col_names=col_names,var_names=var_names, replicate_numbers=replicate_numbers, Conditions=Conditions)) 
}

main_table_read <- function(name_table, extension, working_directory, expression_threshold){

df_final <- data.frame()
df_final_gene <- data.frame()
concat_df <- data.frame()
gene_name_list <- character()
expressed_gene_name_list <- character()

#file
#Load data without header
res <- getting_var_names(extension, working_directory)
col_names <- res$col_names
var_names <- res$var_names

main_table <- data.frame()

main_table <- read.delim(paste0(working_directory,"/",name_table),
                      header = FALSE, # set header = FALSE to avoid overwriting it
                      sep="\t",
                      col.names = col_names )


for (gene in unique(main_table$gene)){
  gene_name_list <- c(gene_name_list, gene)
}
sorted_list <- sort(gene_name_list, decreasing = F)
#print(length(sorted_list))
#####

# Get the column names with the suffix "score"
score_columns <- grep("score$", col_names, value = TRUE)

# Initialize an empty list to store the mean column names
mean_column_names <- list()
expressed_gene_name <- data.frame()
expressed_plus <- data.frame()

expressed_transcript_name <- main_table %>%
  group_by(transcript) %>%
  dplyr::summarize(gene=gene[1],strand=strand[1],
                   across(all_of(score_columns), ~ mean(., na.rm = TRUE), .names = "{.col}_mean")
  )


expressed_plus <- expressed_transcript_name %>%
  filter(strand=="+") %>% 
  select(gene, transcript, strand, contains("plus"))  %>%
  filter(across(all_of(contains("score")), ~ !is.na(.) ) ) %>%
  filter(across(all_of(contains("score")), ~ . >expression_threshold ) )

expressed_minus <- expressed_transcript_name %>%
  filter(strand=="-") %>% 
  select(gene,transcript, strand, contains("minus")) %>%
  filter(across(all_of(contains("score")), ~ !is.na(.) ) ) %>%
  filter(across(all_of(contains("score")), ~ . >expression_threshold ) )

expressed_transcript_name_list <- bind_rows(expressed_plus, expressed_minus) %>% arrange(transcript) %>% pull(transcript)
#print(length(expressed_transcript_name_list)) ## this will set the loading bar

return(list(main_table=main_table,expressed_transcript_name_list=expressed_transcript_name_list))

}

library(tidyr)
library(purrr)
library(dplyr)

results_main_table <- main_table_read(name_table, extension, working_directory, 0.1)
saveRDS(results_main_table, file = "/g/romebioinfo/Projects/tepr/testfromscratch/results_main_table.rds")

message("The results obtained have the following features:")
print(is(results_main_table))
print(length(results_main_table))
print(is(results_main_table[[1]]))
print(is(results_main_table[[2]]))
print(dim(results_main_table[[1]]))
print(length(results_main_table[[2]]))
print(head(results_main_table[[1]]))
print(head(results_main_table[[2]]))
```

The code above was copied to `results_main_table.R` and run with:

```
Rscript results_main_table.R
```

The script should output:

```

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Warning messages:
1: Using `across()` in `filter()` was deprecated in dplyr 1.0.8.
ℹ Please use `if_any()` or `if_all()` instead.
2: Using `across()` in `filter()` was deprecated in dplyr 1.0.8.
ℹ Please use `if_any()` or `if_all()` instead.
3: Using `across()` in `filter()` was deprecated in dplyr 1.0.8.
ℹ Please use `if_any()` or `if_all()` instead.
4: Using `across()` in `filter()` was deprecated in dplyr 1.0.8.
ℹ Please use `if_any()` or `if_all()` instead.
The results obtained have the following features:
[1] "list"   "vector"
[1] 2
[1] "data.frame" "list"       "oldClass"   "vector"
[1] "character"           "vector"              "data.frameRowLabels"
[4] "SuperClassMethod"
[1] 6493200      25
[1] 15002
  biotype  chr  coor1  coor2        transcript   gene strand window
1  lncRNA chr1 817371 817383 ENST00000326734.2 FAM87B      +      1
2  lncRNA chr1 817383 817395 ENST00000326734.2 FAM87B      +      2
3  lncRNA chr1 817395 817407 ENST00000326734.2 FAM87B      +      3
4  lncRNA chr1 817407 817419 ENST00000326734.2 FAM87B      +      4
5  lncRNA chr1 817419 817431 ENST00000326734.2 FAM87B      +      5
6  lncRNA chr1 817431 817443 ENST00000326734.2 FAM87B      +      6
                            id    ctrl_rep1.plus ctrl_rep1.plus_score
1 ENST00000326734.2_FAM87B_+_1 ctrl_rep1.forward                    0
2 ENST00000326734.2_FAM87B_+_2 ctrl_rep1.forward                   NA
3 ENST00000326734.2_FAM87B_+_3 ctrl_rep1.forward                   NA
4 ENST00000326734.2_FAM87B_+_4 ctrl_rep1.forward                   NA
5 ENST00000326734.2_FAM87B_+_5 ctrl_rep1.forward                   NA
6 ENST00000326734.2_FAM87B_+_6 ctrl_rep1.forward                   NA
    ctrl_rep1.minus ctrl_rep1.minus_score    ctrl_rep2.plus
1 ctrl_rep1.reverse                    NA ctrl_rep2.forward
2 ctrl_rep1.reverse                    NA ctrl_rep2.forward
3 ctrl_rep1.reverse                    NA ctrl_rep2.forward
4 ctrl_rep1.reverse                    NA ctrl_rep2.forward
5 ctrl_rep1.reverse                    NA ctrl_rep2.forward
6 ctrl_rep1.reverse                    NA ctrl_rep2.forward
  ctrl_rep2.plus_score   ctrl_rep2.minus ctrl_rep2.minus_score    HS_rep1.plus
1                    0 ctrl_rep2.reverse                    NA HS_rep1.forward
2                   NA ctrl_rep2.reverse                    NA HS_rep1.forward
3                   NA ctrl_rep2.reverse                    NA HS_rep1.forward
4                   NA ctrl_rep2.reverse                    NA HS_rep1.forward
5                   NA ctrl_rep2.reverse                    NA HS_rep1.forward
6                   NA ctrl_rep2.reverse                    NA HS_rep1.forward
  HS_rep1.plus_score   HS_rep1.minus HS_rep1.minus_score    HS_rep2.plus
1                  0 HS_rep1.reverse                  NA HS_rep2.forward
2                 NA HS_rep1.reverse                  NA HS_rep2.forward
3                 NA HS_rep1.reverse                  NA HS_rep2.forward
4                 NA HS_rep1.reverse                  NA HS_rep2.forward
5                 NA HS_rep1.reverse                  NA HS_rep2.forward
6                 NA HS_rep1.reverse                  NA HS_rep2.forward
  HS_rep2.plus_score   HS_rep2.minus HS_rep2.minus_score
1                  0 HS_rep2.reverse                  NA
2                 NA HS_rep2.reverse                  NA
3                 NA HS_rep2.reverse                  NA
4                 NA HS_rep2.reverse                  NA
5                 NA HS_rep2.reverse                  NA
6                 NA HS_rep2.reverse                  NA
[1] "ENST00000000233.10" "ENST00000000412.8"  "ENST00000000442.11"
[4] "ENST00000001008.6"  "ENST00000001146.7"  "ENST00000002125.9"
```

The following code computes ECDF on the main table:

```
results_main_table <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/results_main_table.rds")
main_table <- results_main_table$main_table # %>% filter(gene=="DAP")


```