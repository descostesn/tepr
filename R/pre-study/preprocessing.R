####################
# This script aims at performing pre-processing steps using R only.
#
# Descostes - June 2024 - R-4.4.1
####################


library("AnnotationHub")
library("GenomeInfoDb")
library("GenomicRanges")
library("rtracklayer")
library("parallel")
#library("excluderanges")



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

## Files obtained with bash
protcodbedshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/MANE_Select.protein_coding.bed" # nolint
lncrnabedshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/Ensembl_canonical_TSL123.lncRNA.bed" # nolint
protcodbednoblackwindshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/makewindow/v43.MANE_protein.window200.bed" # nolint
lncrnanednoblackwindshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed" # nolint
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
protcodnoblackfromshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/tmp2.bed" # nolint
lncrnanoblackfromshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/tmp4.bed" # nolint


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
    return(df)
}

makewindowsbedtools <- function(expgr, binsize) {

    ## Filtering out intervals smaller than binsize
    idxsmall <- which(GenomicRanges::width(expgr) < binsize)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        message("Excluding ", lsmall, "/", length(expgr), " annotations that ",
        "are too short.")
        expgr <- expgr[-idxsmall]
    }
    ## command retrieved with HelloRanges:
    ## bedtools_makewindows("-n 200 -b stdin.bed") # nolint
    ## Note: In R, bedtools does not have the "-i srcwinnum" option
    res <- GenomicRanges::tile(expgr, n = binsize)
    res <- unlist(res, use.names = FALSE)

    ## Making names of each element of the list unique
    names(res) <- make.unique(names(res), sep = "_frame")
    return(res)
}

excludeorkeepgrlist <- function(expgr, removegr, removefrom = TRUE,
    ignorestrand = TRUE) {
    ## command retrieved with HelloRanges:
    # nolint - bedtools_intersect("-a protcod.bed -b hg38-blacklist.v2.bed -v")
    resgr <- IRanges::subsetByOverlaps(expgr, removegr,
        invert = removefrom, ignore.strand = ignorestrand)
    return(resgr)
}

bedtogr <- function(currentbed, strand = TRUE) {

    grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
            ranges = IRanges::IRanges(start = currentbed[, 2],
                                  end = currentbed[, 3],
                                  names = currentbed[, 4]),
            strand = if (strand) currentbed[, 6] else "+")
    return(grres)
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




## ------------------------------------------------------------------
## REMOVE


verifybed <- function(bed1, bed2, nbcol = 6) {

    if (!isTRUE(all.equal(nrow(bed1), nrow(bed2))))
        stop("bed1 and bed2 have different nb of rows")

    bedstr1 <- paste(bed1[, 1], bed1[, 2], bed1[, 3], bed1[, 4], bed1[, 5], if(nbcol == 6) bed1[, 6], sep="-") # nolint
    bedstr2 <- paste(bed2[, 1], bed2[, 2], bed2[, 3], bed2[, 4], bed2[, 5], if(nbcol == 6) bed2[, 6], sep="-") # nolint
    idx <- match(bedstr1, bedstr2)
    idxna <- which(is.na(idx))
    lna <- length(idxna)
    if (!isTRUE(all.equal(lna, 0)))
        stop("The gene symbols-chrom are different")
    bed2 <- bed2[idx, ]
    invisible(sapply(seq_len(5), function(i) {
        idx <- which(bed1[, i] != bed2[, i])
        lidx <- length(idx)
        if (!isTRUE(all.equal(lidx, 0)))
            stop("Difference in col ", i)
    }))
}

removepary <- function(dfbed) {
    idx <- grep("_PAR_Y", dfbed[, 4])
    if (!isTRUE(all.equal(length(idx), 0))) {
        message("\t\t Remove PARY")
        dfbed[idx, 4] <- gsub("_PAR_Y", "", dfbed[idx, 4])
    }
    return(dfbed)
}

comparenoblack <- function(bashpath, dfbed) {
    ## Read file obtained with bash
    fromsh <- read.delim(bashpath, header = FALSE)

    ## Remove suffix "_PAR_Y" if present in dfbed
    dfbed <- removepary(dfbed)

    ## Remove last column and add transcript names and strand to match the robj
    ## format
    tmplist <- strsplit(fromsh[, 4], "_")
    tmpnames <- sapply(tmplist, "[", 1)
    tmpstrand <- sapply(tmplist, function(x) x[length(x)])
    fromsh <- data.frame(fromsh$V1, fromsh$V2, fromsh$V3, tmpnames, tmpstrand)

    verifybed(dfbed, fromsh, nbcol = 5)
    return(list(fromsh, dfbed))
}

grtobed <- function(grobj) {
    res <- data.frame(seqnames(grobj), start(grobj), end(grobj), names(grobj),
    strand = strand(grobj))
    return(res)
}

separateframe <- function(dfbed) {

        reslist <- strsplit(dfbed[, 4], "_")
        trsnamevec <- sapply(reslist,  "[", 1)

        if (isTRUE(all.equal(ncol(dfbed), 5))) { # bed from r

            framevec <- sapply(reslist, "[", 2)
            strandvec <- dfbed[, 5]

            ## Makes framevec 1-based and numeric
            framevec <- as.numeric(gsub("frame", "", framevec))
            framevec[which(is.na(framevec))] <- 0
            framevec <- framevec + 1
        } else {
            framevec <- as.numeric(sapply(reslist, function(x) x[length(x)]))
            strandvec <- sapply(reslist, function(x) x[length(x) - 1])
        }

        ## Verify there is one frame per transcript line
        if (!isTRUE(all.equal(length(trsnamevec), length(framevec))))
            stop("Problem of correspondance between trsnamevec and framevec")

        ## Separate transcript names from frame in two columns 
        resbed <- data.frame(chrom = dfbed[, 1], start = dfbed[, 2],
            end = dfbed[, 3], name = trsnamevec, strand = strandvec,
            frame = framevec)
        return(resbed)
}

comparewind <- function(fromr_noblackshgr, fromsh_noblackwindpath, windsize) {
    ## Preparing bed df
    fomr_windgr <- makewindowsbedtools(fromr_noblackshgr, windsize)
    fromr_windbed <- grtobed(fomr_windgr)
    fromsh_windbed <- read.delim(fromsh_noblackwindpath, header = FALSE)

    ## Remove suffix "_PAR_Y" if present in bed
    fromr_windbed <- removepary(fromr_windbed)
    fromsh_windbed <- removepary(fromsh_windbed)

    ## Separate transcript names from frame in two columns
    fromr_windbed <- separateframe(fromr_windbed)
    fromsh_windbed <- separateframe(fromsh_windbed)

    verifybed(fromr_windbed, fromsh_windbed)
}

## Compare the bed files before removing black lists
protcodbedsh <- read.delim(protcodbedshpath, header = FALSE)
lncrnabedsh <- read.delim(lncrnabedshpath, header = FALSE)
verifybed(protcodbed, protcodbedsh)
verifybed(lncrnabed, lncrnabedsh)

## Exclude black list with the file that was used in bash
blacklistsh <- read.delim(blacklistshpath, header = FALSE)
blacklistshgr <- bedtogr(blacklistsh, strand = FALSE)
protcodnoblackshgr <- excludeorkeepgrlist(protcodgr, blacklistshgr)
lncrnanoblackshgr <- excludeorkeepgrlist(lncrnagr, blacklistshgr)
protcodnoblacksh <- grtobed(protcodnoblackshgr)
lncrnanoblacksh <- grtobed(lncrnanoblackshgr)

## Compare protcodnoblacksh and lncrnanoblacksh to the files obtained with
## bash
res <- comparenoblack(protcodnoblackfromshpath, protcodnoblacksh)
fromsh_protcodnoblackshbed <- res[[1]]
fromr_protcodnoblackshbed <- res[[2]]
res <- comparenoblack(lncrnanoblackfromshpath, lncrnanoblacksh)
fromsh_lncrnanoblackshbed <- res[[1]]
fromr_lncrnanoblackshbed <- res[[2]]

## Temporary variables for comparison with files obtained with bash
fomr_protcodwindgr <- makewindowsbedtools(protcodnoblackshgr, windsize)
fromr_lncrnawindgr <- makewindowsbedtools(lncrnanoblackshgr, windsize)
fomr_protcodwindbed <- grtobed(fomr_protcodwindgr)
fromr_lncrnawindbed <- grtobed(fromr_lncrnawindgr)

## Compare with bash files
fromsh_protcodwindbed <- read.delim(protcodbednoblackwindshpath, header = FALSE)
fromsh_lncrnawindbed <- read.delim(lncrnanednoblackwindshpath, header = FALSE)
comparewind(protcodnoblackshgr, protcodbednoblackwindshpath, windsize)
comparewind(lncrnanoblackshgr, lncrnanednoblackwindshpath, windsize)


################## DECIPHERING  CODE
# buildstr <- function(gencode) {
#     trsidvec <-unlist(lapply(strsplit(unlist(lapply(strsplit(gencode$V9, ";"),"[",2))," "),"[",3))
# symbolvec <-unlist(lapply(strsplit(unlist(lapply(strsplit(gencode$V9, ";"),"[",4))," "),"[",3))
# gencodestr <- paste(gencode$V1, gencode$V4, gencode$V5, trsidvec, symbolvec, gencode$V7,sep="-")
# return(gencodestr)
# }
# tofind <- "chrX-104072887-104076236-ENST00000598087.4-TMSB15B--"
# strcomp <- buildstr(gencode)
# grep(tofind, strcomp)
# gencodestrman <- read.delim("/g/romebioinfo/Projects/tepr/downloads/annotations/temp1.gtf", header=F)
# strcomp <- buildstr(gencodestrman)
# grep(tofind, strcomp)
# gencodestrman <- read.delim("/g/romebioinfo/Projects/tepr/downloads/annotations/tmp2.bed", header=F)
# strcomp <- paste(gencodestrman$V1, gencodestrman$V2, gencodestrman$V3, gsub("_","-", gencodestrman$V4), sep="-")
# grep(tofind, strcomp)
# protcodbedsh <- read.delim(protcodbedshpath, header = FALSE)
# strcomp <- paste(protcodbedsh$V1, protcodbedsh$V2, protcodbedsh$V3, protcodbedsh$V4, protcodbedsh$V5, protcodbedsh$V6, sep="-")
# grep(tofind, strcomp)
# strcomp <- buildstr(gencodeprotcod)
# grep(tofind, strcomp)
# strcomp <- paste(protcodbed$chrom, protcodbed$start, protcodbed$end, protcodbed$ensembl, protcodbed$symbol, protcodbed$strand, sep="-")
# grep(tofind, strcomp)
# !!!!!!!!!!!!!!!!!!!!!!!
# comparenoblack(protcodnoblackfromshpath, protcodnoblacksh)
# ## Check if the annotations were present before making the windows in bash
# protcodnoblackfromsh <- read.delim(protcodnoblackfromshpath, header = FALSE)
# annotrsnoblackfromsh <- unique(sapply(strsplit(protcodnoblackfromsh$V4, "_"), "[", 1))
# !!!!!!!!!!!!!!!!!!!!!!!!
# ## Verify transcripts annotations between r and sh
# annotrsfromsh <- unique(sapply(strsplit(protcodwindfromsh$V4, "_"), "[", 1))
# annotrsfromr <- unique(sapply(strsplit(names(protcodwindowstmp), "_"), "[", 1))
# idx <- match(annotrsfromr, annotrsfromsh)
# ## The r variable carries more annotations
# idxna <- which(is.na(idx))
# head(annotrsfromr[idxna])


## End REMOVE
## ------------------------------------------------------------------




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
