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
library("clusterProfiler")
library("org.Hs.eg.db")

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
nbcpu <- 6
database_name <- "org.Hs.eg.db"


##################
#FUNCTIONS
##################

.returnwindowvec <- function(dfintervalsrownames) {
        windowvec <- as.numeric(
        gsub("frame", "",
            sapply(
                strsplit(dfintervalsrownames, "_"),
            "[", 2)))
    windowvec[which(is.na(windowvec))] <- 0
    windowvec <- windowvec + 1
    return(windowvec)
}

.returnsymbolvec <- function(transvec, database_name, dfintervals) {

    transrootvec <- gsub("\\..+", "", transvec, perl = TRUE)
    transidvec <- unique(transrootvec)

    convertdf <- clusterProfiler::bitr(transidvec, fromType = "ENSEMBLTRANS", 
            toType = c("ENSEMBLTRANS", "SYMBOL"), OrgDb = database_name)
    idxdup <- which(duplicated(convertdf[,1]))
    if (!isTRUE(all.equal(length(idxdup), 0)))
        convertdf <- convertdf[-idxdup, ]

    #### Corresponding values
    idx <- match(transrootvec, convertdf[, 1])
    symbolvec <- convertdf[idx, 2]
    if (!isTRUE(all.equal(length(symbolvec), nrow(dfintervals))))
        stop("Problem in retrieving symbols")

    return(symbolvec)
}

.retrievemeanfrombw <- function(grintervals, bwpath, verbose) {

    rangeselect <- rtracklayer::BigWigSelection(grintervals, character())
    bwval <- rtracklayer::import.bw(bwpath,
        selection = rangeselect, as = "NumericList")

    if (!isTRUE(all.equal(length(bwval), length(grintervals))))
        stop("The number of intervals retrieved from the bigwig is not correct")
    if (!isTRUE(all.equal(names(bwval), names(grintervals)))) {
        if (verbose) message("\t Re-ordering list")
        idx <- match(names(grintervals), names(bwval))
        idxna <- which(is.na(idx))
        lna <- length(idxna)
        if (!isTRUE(all.equal(lna, 0)))
            stop("Problem with matching names.")
        bwval <- bwval[idx]
    }

    if (verbose) message("\t Computing mean values for each interval")
    meanvec <- sapply(bwval, mean)
    rm(bwval)
    invisible(gc())
    return(meanvec)
}

buildscoreforintervals <- function(grintervals, expdf, grname, nbcpu,
    database_name, verbose = TRUE) {

    if (verbose) message("Processing ", grname)

    scorelist <- mcmapply(function(bwpath, expname, grintervals, verbose) {
        if (verbose) message("\t Retrieving bw values for ", expname)
        meanvec <- .retrievemeanfrombw(grintervals, bwpath, verbose)
        return(meanvec)
    }, expdf$path, expdf$name, MoreArgs = list(grintervals, verbose),
        SIMPLIFY = FALSE, mc.cores = nbcpu)

    ## Building matrix of mean scores
    if (verbose) message("\t Building matrix")
    scoremat <- do.call("cbind", scorelist)
    colnames(scoremat) <- expdf$name
    dfintervals <- as.data.frame(grintervals)
    if (!isTRUE(all.equal(nrow(scoremat), nrow(dfintervals))))
        stop("Differing number of rows for score matrix and annotations in ",
            "function buildscoreforintervals")
    if (!isTRUE(all.equal(rownames(scoremat), rownames(dfintervals))))
        stop("The rows of scoremat and dfintervals are not in the same order")

    ## Building vectors used for the final data.frame
    if (verbose) message("\t Retrieving information for the final data.frame")
    dfintervalsrownames <- rownames(dfintervals)
    transvec <- gsub("_frame.+", "", dfintervalsrownames, perl = TRUE) # nolint
    symbolvec <- .returnsymbolvec(transvec, database_name, dfintervals)
    windowvec <- .returnwindowvec(dfintervalsrownames)
    strandvec <- dfintervals$strand

    ## Final data.frame
    if (verbose) message("\t Creating the final data.frame")
    df <- data.frame(biotype = grname, chr = dfintervals$seqnames,
    start = dfintervals$start, end = dfintervals$end, transcript = transvec,
    gene = symbolvec, strand = strandvec, window = windowvec,
    id = paste(transvec, symbolvec, strandvec, windowvec, sep = "_"))
    df <- cbind(df, scoremat)

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
    nbcpu, database_name)
lncrnadf <- buildscoreforintervals(lncrnawindows, exptab, "lncrna", nbcpu,
    database_name)
alldf <- rbind(protcoddf, lncrnadf)

