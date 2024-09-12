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


