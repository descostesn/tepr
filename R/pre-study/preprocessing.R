####################
# This script aims at performing the pre-processing steps using R only.
#
# Descostes - June 2024 - R-4.4.1
####################


library("AnnotationHub")
library("GenomeInfoDb")
library("GenomicRanges")
library("rtracklayer")
library("parallel")

source("commons.R")


##################
# PARAMETERS
##################

gencodepath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
## Note: For a complete list of blacklist names see
## ah <- AnnotationHub() # nolint
## query_data <- subset(ah, preparerclass == "excluderanges") # nolint
## print(query_data) # nolint
blacklistname <- "hg38.Kundaje.GRCh38_unified_Excludable"
outputfolder <- "/g/romebioinfo/Projects/tepr/downloads"
robjoutputfold <- "/g/romebioinfo/Projects/tepr/robjsave"
## Mappability tracks in bed format can be downloaded from https://bismap.hoffmanlab.org/ # nolint
## Scroll down to the table containing pre-computed tracks for hg38, hg19, mm10,
## and mm9 in bins of 24, 36, 50, and 100 bp in single- or multireads.
## The bed file below is for hg38 unique reads of 50 bp.
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.Unique.Mappability.bed" # nolint
## Size of the window to extract values
windsize <- 200
## Table of experiments - contains the columns "name,condition,replicate,strand,path" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv"
nbcpu <- 15



##################
#FUNCTIONS
##################

retrievemeanfrombw <- function(grintervals, bwpath) {
    message("Retrieving bw values")
    rangeselect <- rtracklayer::BigWigSelection(grintervals, character())
    bwval <- rtracklayer::import.bw(bwpath,
        selection = rangeselect, as = "NumericList")

    if (!isTRUE(all.equal(length(bwval), length(grintervals))))
        stop("The number of intervals retrieved from the bigwig is not correct")
    if (!isTRUE(all.equal(names(bwval), names(grintervals)))) {
        message("\t Re-ordering list")
        idx <- match(names(grintervals), names(bwval))
        idxna <- which(is.na(idx))
        lna <- length(idxna)
        if (!isTRUE(all.equal(lna, 0)))
            stop("Problem with matching names.")
        bwval <- bwval[idx]
    }

    message("\t Computing mean values for each interval")
    meanvec <- sapply(bwval, mean)
    rm(bwval)
    invisible(gc())
    return(meanvec)
}

buildscoreforintervals <- function(grintervals, expdf, grname, nbcpu) {

    message("Retrieving values for ", grname)

    scorelist <- mcmapply(function(bwpath, expname, grintervals) {
        message("\t Retrieving values for ", expname)
        meanvec <- retrievemeanfrombw(grintervals, bwpath)
        return(meanvec)
    }, expdf$path, expdf$name, MoreArgs = list(grintervals), SIMPLIFY = FALSE,
        mc.cores = nbcpu)

    scoremat <- do.call("cbind", scorelist)
    colnames(scoremat) <- expdf$name
    dfintervals <- as.data.frame(grintervals)
    if (!isTRUE(all.equal(nrow(scoremat), nrow(dfintervals))))
        stop("Differing number of rows for score matrix and annotations in ",
            "function buildscoreforintervals")
    if (!isTRUE(all.equal(rownames(scoremat), rownames(dfintervals))))
        stop("The rows of scoremat and dfintervals are not in the same order")
    df <- cbind(dfintervals, scoremat)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    transvec <- gsub("_frame.+", "", rownames(dfintervals), perl = TRUE) # nolint
    df <- data.frame(biotype = grname, chr = dfintervals$seqnames,
        start = dfintervals$start, end = dfintervals$end, transcript = transvec,
        gene = 
biotype  chr  coor1  coor2        transcript   gene strand window
1 protein-coding chr1 923923 924026 ENST00000616016.5 SAMD11      +      1
2 protein-coding chr1 924026 924129 ENST00000616016.5 SAMD11      +      2
3 protein-coding chr1 924129 924232 ENST00000616016.5 SAMD11      +      3
4 protein-coding chr1 924232 924335 ENST00000616016.5 SAMD11      +      4
5 protein-coding chr1 924335 924438 ENST00000616016.5 SAMD11      +      5
6 protein-coding chr1 924438 924541 ENST00000616016.5 SAMD11      +      6
                            id         dataset.x  score.x         dataset.y
1 ENST00000616016.5_SAMD11_+_1 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
2 ENST00000616016.5_SAMD11_+_2 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
3 ENST00000616016.5_SAMD11_+_3 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
4 ENST00000616016.5_SAMD11_+_4 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
5 ENST00000616016.5_SAMD11_+_5 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
6 ENST00000616016.5_SAMD11_+_6 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
  score.y       dataset.x.x score.x.x       dataset.y.y score.y.y
1    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
2    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
3    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
4    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
5    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
6    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
    dataset.x.x.x score.x.x.x   dataset.y.y.y score.y.y.y dataset.x.x.x.x
1 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
2 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
3 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
4 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
5 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
6 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
  score.x.x.x.x dataset.y.y.y.y score.y.y.y.y
1      0.000000 HS_rep2.reverse          <NA>
2      0.000000 HS_rep2.reverse          <NA>
3      0.000000 HS_rep2.reverse          <NA>
4      0.000000 HS_rep2.reverse          <NA>
5      0.000000 HS_rep2.reverse          <NA>
6      0.000000 HS_rep2.reverse          <NA>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    return(df)
}


createfolder <- function(outfold) {
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
}

createblacklist <- function(blacklistname, outputfolder) { # nolint

    blacklistgr <- AnnotationHub::query(AnnotationHub::AnnotationHub(),
        blacklistname)[[1]]
    blacklistgr <- blacklistgr %>%
                   sort() %>%
                   GenomeInfoDb::keepStandardChromosomes(pruning.mode = "tidy")
    createfolder(outputfolder)
    write.table(as.data.frame(blacklistgr),
            file = file.path(outputfolder, paste0(blacklistname, ".bed")),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(blacklistgr)
}

grepsequential <- function(valvec, gentab, invert = FALSE, verbose = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ]
    }))
    return(gentab)
}

sortedbedformat <- function(gencode) {
    gencode <- gencode[order(gencode$V1, gencode$V4), ] # nolint ## Ordering by chrom and start
    infolist <- strsplit(gencode$V9, ";")
    namevec <- gsub(" gene_name ", "", sapply(infolist, "[", 4)) # nolint
    ensnamevec <- gsub(" transcript_id ", "", sapply(infolist, "[", 2)) # nolint
    gencodebed <- cbind(gencode[, c(1, 4, 5)], ensnamevec, namevec,
        gencode[, 7])
    colnames(gencodebed) <- c("chrom", "start", "end", "ensembl", "symbol",
        "strand")
    return(gencodebed)
}


##################
# MAIN
##################

createfolder(robjoutputfold)

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)
gencode <- gencode[which(gencode$V3 == "transcript"), ]

## Selecting Ensembl_canonical transcripts i.e. most representative transcript
## of the protein coding gene. This will be the MANE_Select transcript if there
## is one, or a transcript chosen by an Ensembl algorithm otherwise.
gencodeprotcod <- grepsequential("MANE_Select", gencode)
protcodbed <- sortedbedformat(gencodeprotcod)
protcodgr <- bedtogr(protcodbed)

## Retrieve long non-coding transcripts
lncrna <- grepsequential(c("lncRNA", "Ensembl_canonical"), gencode)
removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
lncrna <- grepsequential(removevec, lncrna, invert = TRUE)
lncrnabed <- sortedbedformat(lncrna)
lncrnagr <- bedtogr(lncrnabed)

## Exclude blacklist
blacklistgr <- createblacklist(blacklistname, outputfolder)
protcodnoblackgr <- excludeorkeepgrlist(protcodgr, blacklistgr)
lncrnanoblackgr <- excludeorkeepgrlist(lncrnagr, blacklistgr)

## Saving objects to check conformity with bash results
saveRDS(protcodbed, file = file.path(robjoutputfold, "protcodbed.rds"))
saveRDS(protcodgr, file = file.path(robjoutputfold, "protcodgr.rds"))
saveRDS(lncrnabed, file = file.path(robjoutputfold, "lncrnabed.rds"))
saveRDS(lncrnagr, file = file.path(robjoutputfold, "lncrnagr.rds"))

## Exclude low mappability
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES - the mappability track
## used has only 1 as mapping scores. See parameters.
maptrack <- read.delim(maptrackpath, header = FALSE)
maptrackgr <- bedtogr(maptrack)
protcodnoblacknomapgr <- excludeorkeepgrlist(protcodnoblackgr, maptrackgr,
    removefrom = FALSE)
lncrnanoblacknomapgr <- excludeorkeepgrlist(lncrnanoblackgr, maptrackgr,
    removefrom = FALSE)

## Make windows of windsize for each annotation
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES
protcodwindows <- makewindowsbedtools(protcodnoblacknomapgr, windsize)
lncrnawindows <- makewindowsbedtools(lncrnanoblacknomapgr, windsize)

## Retrieving values from bigwig files
exptab <- read.csv(exptabpath, header = TRUE)
protcoddf <- buildscoreforintervals(protcodwindows, exptab, "protein_coding",
    nbcpu)
lncrnadf <- buildscoreforintervals(lncrnawindows, exptab, "lncrna", nbcpu)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
check in above variables the columns compared to what is below and then merge
the two df. Save the object and compare it to the df_bound object

!! TODO: the last filter remove the PAR genes ie.e pseudoautosomal genes both
!! in X and Y by doing 'strand != "Y"'

!!The data frame should contain the columns
[1] "biotype"               "chr"                   "coor1"
 [4] "coor2"                 "transcript"            "gene"
 [7] "strand"                "window"                "id"
[10] "ctrl_rep1.plus"        "ctrl_rep1.plus_score"  "ctrl_rep1.minus"
[13] "ctrl_rep1.minus_score" "ctrl_rep2.plus"        "ctrl_rep2.plus_score"
[16] "ctrl_rep2.minus"       "ctrl_rep2.minus_score" "HS_rep1.plus"
[19] "HS_rep1.plus_score"    "HS_rep1.minus"         "HS_rep1.minus_score"
[22] "HS_rep2.plus"          "HS_rep2.plus_score"    "HS_rep2.minus"
[25] "HS_rep2.minus_score"

Here is an example:

[[8]]
         biotype  chr  coor1  coor2        transcript   gene strand window
1 protein-coding chr1 923923 924026 ENST00000616016.5 SAMD11      +      1
2 protein-coding chr1 924026 924129 ENST00000616016.5 SAMD11      +      2
3 protein-coding chr1 924129 924232 ENST00000616016.5 SAMD11      +      3
4 protein-coding chr1 924232 924335 ENST00000616016.5 SAMD11      +      4
5 protein-coding chr1 924335 924438 ENST00000616016.5 SAMD11      +      5
6 protein-coding chr1 924438 924541 ENST00000616016.5 SAMD11      +      6
                            id         dataset score
1 ENST00000616016.5_SAMD11_+_1 HS_rep2.reverse  <NA>
2 ENST00000616016.5_SAMD11_+_2 HS_rep2.reverse  <NA>
3 ENST00000616016.5_SAMD11_+_3 HS_rep2.reverse  <NA>
4 ENST00000616016.5_SAMD11_+_4 HS_rep2.reverse  <NA>
5 ENST00000616016.5_SAMD11_+_5 HS_rep2.reverse  <NA>
6 ENST00000616016.5_SAMD11_+_6 HS_rep2.reverse  <NA>

> lapply(list_of_dfs, function(x) unique(x$biotype))
[[1]]
[1] "protein-coding"

[[2]]
[1] "protein-coding"

[[3]]
[1] "protein-coding"

[[4]]
[1] "protein-coding"

[[5]]
[1] "protein-coding"

[[6]]
[1] "protein-coding"

[[7]]
[1] "protein-coding"

[[8]]
[1] "protein-coding"

After retrieving scores for all files we obtain the variable joined_df (see
joined_df path in parameters of R\pre-study\conformity_with_bash.R) with:

 biotype  chr  coor1  coor2        transcript   gene strand window
1 protein-coding chr1 923923 924026 ENST00000616016.5 SAMD11      +      1
2 protein-coding chr1 924026 924129 ENST00000616016.5 SAMD11      +      2
3 protein-coding chr1 924129 924232 ENST00000616016.5 SAMD11      +      3
4 protein-coding chr1 924232 924335 ENST00000616016.5 SAMD11      +      4
5 protein-coding chr1 924335 924438 ENST00000616016.5 SAMD11      +      5
6 protein-coding chr1 924438 924541 ENST00000616016.5 SAMD11      +      6
                            id         dataset.x  score.x         dataset.y
1 ENST00000616016.5_SAMD11_+_1 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
2 ENST00000616016.5_SAMD11_+_2 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
3 ENST00000616016.5_SAMD11_+_3 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
4 ENST00000616016.5_SAMD11_+_4 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
5 ENST00000616016.5_SAMD11_+_5 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
6 ENST00000616016.5_SAMD11_+_6 ctrl_rep1.forward 0.000000 ctrl_rep1.reverse
  score.y       dataset.x.x score.x.x       dataset.y.y score.y.y
1    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
2    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
3    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
4    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
5    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
6    <NA> ctrl_rep2.forward  0.000000 ctrl_rep2.reverse      <NA>
    dataset.x.x.x score.x.x.x   dataset.y.y.y score.y.y.y dataset.x.x.x.x
1 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
2 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
3 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
4 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
5 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
6 HS_rep1.forward    0.000000 HS_rep1.reverse        <NA> HS_rep2.forward
  score.x.x.x.x dataset.y.y.y.y score.y.y.y.y
1      0.000000 HS_rep2.reverse          <NA>
2      0.000000 HS_rep2.reverse          <NA>
3      0.000000 HS_rep2.reverse          <NA>
4      0.000000 HS_rep2.reverse          <NA>
5      0.000000 HS_rep2.reverse          <NA>
6      0.000000 HS_rep2.reverse          <NA>

