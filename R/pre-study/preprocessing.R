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

# source("commons.R")


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
## The bed file below was created and sent by Victor
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
## Size of the window to extract values
windsize <- 200
## Table of experiments - contains the columns "name,condition,replicate,strand,path" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv"
nbcpu <- 6
database_name <- "org.Hs.eg.db"


##################
#FUNCTIONS
##################


createfolder <- function(outfold) {
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
}

checkexptab <- function(exptab) {
    colnamevec <- c("condition", "direction", "path", "replicate", "strand")
    if (!isTRUE(all.equal(sort(colnames(exptab)), colnamevec)))
        stop("The experiment table should have the columns: ",
            "'condition', 'direction', 'path', 'replicate', 'strand'")

    if (!isTRUE(all.equal(length(unique(exptab$condition)), 2)))
        stop("The table should only contain two conditions")

    if (isTRUE(all.equal(length(grep("ctrl", exptab$condition)), 0)))
        stop("The control condition (or condition 1) should be designated",
            " by 'ctrl'")

    directionvec <- unique(exptab$direction)
    if (!isTRUE(all.equal(length(directionvec), 2)) ||
        !isTRUE(all.equal(directionvec, c("fwd", "rev"))))
        stop("Only two values are allowed for the column direction of the",
            "experiment table, 'fwd' and 'rev'")

    strandvec <- unique(exptab$strand)
    if (!isTRUE(all.equal(strandvec, c("+", "-"))))
        stop("The strand column of the experiment table should only contain",
            " '+' and '-'.")
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

checkremoval <- function(datagr, dataremovedgr, dataname, removename,
    toremovegr, removeopt) {

    message("\n\n Checking operations for ", dataname, " with data of ",
        removename)

    ## Calculating number of elements overlapping toremovegr
    res <- GenomicRanges::findOverlaps(datagr, toremovegr)
    nboverdata <- length(unique(S4Vectors::queryHits(res)))

    if (removeopt) {
        message("The number of elements that should be removed is: ",
            nboverdata)
        subres <- length(datagr) - nboverdata
        message("The number of elements of the resulting object after ",
            "subtraction should be: ", length(datagr), "-", nboverdata, "=",
            subres)
        message("The number of elements in the resulting object is: ",
            length(dataremovedgr))
    } else {
        message("The number of elements of the data to keep is: ",
            nboverdata)
        message("The number of elements before the overlap is: ",
            length(datagr))
        message("The number of elements in the resulting object is: ",
            length(dataremovedgr))
    }
}


bedtogr <- function(currentbed, strand = TRUE, symbol = TRUE) {

    grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
            ranges = IRanges::IRanges(start = currentbed[, 2],
                                  end = currentbed[, 3],
                                  names = currentbed[, 4]),
            strand = if (strand) currentbed[, 6] else "*",
            symbol = if (symbol) currentbed[, 5] else NA)
    return(grres)
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

    ## Change row names to keep the gene symbols
    names(expgr) <- paste(names(expgr), expgr$symbol, sep = "_")

    ## command retrieved with HelloRanges:
    ## bedtools_makewindows("-n 200 -b stdin.bed") # nolint
    ## Note: In R, bedtools does not have the "-i srcwinnum" option
    res <- GenomicRanges::tile(expgr, n = binsize)
    res <- unlist(res, use.names = FALSE)

    ## Adding back metadata from names
    tmplist <- strsplit(names(res), "_")
    transvec <- sapply(tmplist, "[", 1)
    symbolvec <- sapply(tmplist, "[", 2)
    names(res) <- transvec
    S4Vectors::elementMetadata(res)[, "symbol"] <- symbolvec

    ## Making names of each element of the list unique
    names(res) <- make.unique(names(res), sep = "_frame")
    return(res)
}


retrieveandfilterfrombg <- function(exptab, blacklistgr, maptrackgr, nbcpu,
    verbose = TRUE) {

    expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)

    ## Looping on each experiment bw file
    # currentpath <- exptab$path[1]
    # currentname <- expnamevec[1]
    bedgraphgrlist <- mcmapply(function(currentpath, currentname, blacklistgr,
        maptrackgr, nbcpu, verbose) {

        if (verbose) message("\t Retrieving values for ", currentname)
        valgr <- rtracklayer::import.bedGraph(currentpath)

        if (verbose) message("\t\t Filtering out scores in black list ranges")
        resblack <- GenomicRanges::findOverlaps(valgr, blacklistgr,
            ignore.strand = TRUE)
        idxblack <- unique(S4Vectors::queryHits(resblack))
        BiocGenerics::score(valgr)[idxblack] <- NA

        if (verbose) message("\t\t Keeping high mappability scores")
        reshigh <- GenomicRanges::findOverlaps(valgr, maptrackgr,
            ignore.strand = TRUE)
        idxhigh <- unique(S4Vectors::queryHits(reshigh))
        if (isTRUE(all.equal(length(idxhigh), length(valgr))))
            message("Only highly mappable element were found")
        else
            ## Setting the scores of the ranges NOT in idxhigh to NA
            BiocGenerics::score(valgr)[-idxhigh] <- NA

        return(valgr)

    }, exptab$path, expnamevec, MoreArgs = list(blacklistgr,
        maptrackgr, nbcpu, verbose), mc.cores = nbcpu, SIMPLIFY = FALSE)

    return(bedgraphgrlist)
}


##################
# MAIN
##################

createfolder(robjoutputfold)

## Reading the information about experiments
exptab <- read.csv(exptabpath, header = TRUE)
checkexptab(exptab)

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)
gencode <- gencode[which(gencode$V3 == "transcript"), ]

## Selecting Ensembl_canonical transcripts i.e. most representative transcript
## of the protein coding gene. This will be the MANE_Select transcript if there
## is one, or a transcript chosen by an Ensembl algorithm otherwise.
gencodeprotcod <- grepsequential("MANE_Select", gencode)
protcodbed <- sortedbedformat(gencodeprotcod)

## Retrieve long non-coding transcripts
lncrna <- grepsequential(c("lncRNA", "Ensembl_canonical"), gencode)
removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
lncrna <- grepsequential(removevec, lncrna, invert = TRUE)
lncrnabed <- sortedbedformat(lncrna)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (verbose) message("Combine the annotations")
allannobed <- rbind(protcodbed, lncrnabed)
allannogr <- bedtogr(allannobed)

if (verbose) message("Make windows for all annotations")
allwindows <- makewindowsbedtools(allannogr, windsize)

if (verbose) message("Reading the black list")
blacklistgr <- createblacklist(blacklistname, outputfolder)

if (verbose) message("Reading the highly mappable ranges")
maptrack <- read.delim(maptrackpath, header = FALSE)
maptrackgr <- bedtogr(maptrack, strand = FALSE)

## Retrieving the values of the bedgraph files, removing black lists and keeping
## high mappability scores
message("Reading and filtering bedgraphs")
bedgraphgrlist <- retrieveandfilterfrombg(allwindows, exptab, blacklistgr,
    maptrackgr, nbcpu)

## Retrieving values according to annotations and calculate an arithmetic
## weighted mean



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Saving objects to check conformity with bash results
saveRDS(protcodbed, file = file.path(robjoutputfold, "protcodbed.rds"))
saveRDS(protcodgr, file = file.path(robjoutputfold, "protcodgr.rds"))
saveRDS(lncrnabed, file = file.path(robjoutputfold, "lncrnabed.rds"))
saveRDS(lncrnagr, file = file.path(robjoutputfold, "lncrnagr.rds"))
# protcodbed <- readRDS(file.path(robjoutputfold, "protcodbed.rds"))
# protcodgr <- readRDS(file.path(robjoutputfold, "protcodgr.rds"))
# lncrnabed <- readRDS(file.path(robjoutputfold, "lncrnabed.rds"))
# lncrnagr <- readRDS(file.path(robjoutputfold, "lncrnagr.rds"))


## Exclude blacklist
blacklistgr <- createblacklist(blacklistname, outputfolder)
protcodnoblackgr <- excludeorkeepgrlist(protcodgr, blacklistgr)
lncrnanoblackgr <- excludeorkeepgrlist(lncrnagr, blacklistgr)

## Check excluded intervals using blacklist
checkremoval(protcodgr, protcodnoblackgr, "proteincoding", "blacklist",
    blacklistgr, removeopt = TRUE)
checkremoval(lncrnagr, lncrnanoblackgr, "lncrna", "blacklist",
    blacklistgr, removeopt = TRUE)

## Exclude low mappability
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES - the mappability track
## used has only 1 as mapping scores. See parameters.
maptrack <- read.delim(maptrackpath, header = FALSE)
maptrackgr <- bedtogr(maptrack)


!!!!!!!!
keepgrlist <- function(datagr, maptrackgr) {
    res <- GenomicRanges::findOverlaps(datagr, maptrackgr)
    idxtokeep <- unique(S4Vectors::queryHits(res))
    datagr <- datagr[idxtokeep, ]
    return(datagr)
}
protcodnoblackmaponlygr <- keepgrlist(protcodnoblackgr, maptrackgr)
lncrnanoblackmaponlygr <- keepgrlist(lncrnanoblackgr, maptrackgr)

!!!!!!!!!!!

protcodnoblacknomapgr <- excludeorkeepgrlist(protcodnoblackgr, maptrackgr,
    removefrom = FALSE)
lncrnanoblacknomapgr <- excludeorkeepgrlist(lncrnanoblackgr, maptrackgr,
    removefrom = FALSE)

## Check excluded intervals because of low mappability for protein coding
checkremoval(protcodnoblackgr, protcodnoblacknomapgr, "proteincoding",
    "maptrack", maptrackgr, removeopt = FALSE)
checkremoval(lncrnanoblackgr, lncrnanoblacknomapgr, "lncrna", "maptrack",
    maptrackgr, removeopt = FALSE)

## Make windows of windsize for each annotation
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES

lncrnawindows <- makewindowsbedtools(lncrnanoblacknomapgr, windsize)

## Retrieving values from bigwig files
protcoddf <- buildscoreforintervals(protcodwindows, exptab, "protein_coding",
    nbcpu, database_name)
lncrnadf <- buildscoreforintervals(lncrnawindows, exptab, "lncrna", nbcpu,
    database_name)
alldf <- rbind(protcoddf, lncrnadf)
saveRDS(alldf, file = file.path(robjoutputfold, "alldffrompreprocessing.rds"))
#alldf <- readRDS("/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds") # nolint
