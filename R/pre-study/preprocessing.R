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
library("purrr")
library("dplyr")



##################
# PARAMETERS
##################

gencodepath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
windsize <- 200
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv" # nolint
database_name <- "org.Hs.eg.db"
outputfolder <- "/g/romebioinfo/tmp/preprocessing"
robjoutputfold <- outputfolder

## Set this variable to NULL if the online retrieval should be performed
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
## The bed file below was created and sent by Victor
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
## Parallelization on bedgraph files. The maximum should be equal to the number of bedgraph files.  # nolint
nbcpubg <- 8
## Parallelization on transcripts. The maximum should be limited to the capacity of your machine.  # nolint
nbcputrans <- 20

## Note: For a complete list of blacklist names see
## ah <- AnnotationHub() # nolint
## query_data <- subset(ah, preparerclass == "excluderanges") # nolint
## print(query_data) # nolint
blacklistname <- "hg38.Kundaje.GRCh38_unified_Excludable"



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

## For different strings provided in the vector "valvec", perform the filtering
## of gentab on each string using the result of the previous filtering
grepsequential <- function(valvec, gentab, invert = FALSE, verbose = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ] ## This line enables sequential grep
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


.divideannoinwindows <- function(expbed, windcoordvec, nbwindows, nbcputrans) {

    cl <- parallel::makeCluster(nbcputrans)
    windflist <- parallel::parLapply(cl, seq_len(nrow(expbed)),
        function(i, expbed, windcoordvec, nbwindows) {

            currentanno <- expbed[i, ]
            ## Retrieve the necessary gene information
            currentstart <- currentanno$start
            currentend <- currentanno$end
            currentstrand <- currentanno$strand
            windowvec <- windcoordvec

            ## Compute the vector with the size of each window
            lgene <- currentend - currentstart
            windowsize <- round(lgene / nbwindows)
            missingbp <- lgene %% nbwindows
            windsizevec <- rep(windowsize, nbwindows)
            ## Add the missing nb of bp (that is ignore by tile) in the last
            ## element of windsizevec
            if (!isTRUE(all.equal(missingbp, 0)))
                windsizevec[nbwindows] <- windsizevec[nbwindows] + missingbp

            ## Building the start and end vectors using the cummulative sum
            cumsumvec <- cumsum(c(currentstart, windsizevec))
            startvec <- cumsumvec[-length(cumsumvec)]
            endvec <- cumsumvec[-1]
            if (!isTRUE(all.equal(endvec - startvec, windsizevec)))
                stop("Problem in the calculation of windows")

            ## Inverting start, end, and window vectors if strand is negative
            if (isTRUE(all.equal(currentstrand, "-"))) {
                startvec <- rev(startvec)
                endvec <- rev(endvec)
                windowvec <- rev(windcoordvec)
            }

            ## Build the result data.frame containing the coordinates of each
            ## frame alongside window and coord numbers
            res <- data.frame(biotype = currentanno$biotype,
                chr = currentanno$chrom, coor1 = startvec,
                coor2 = endvec,  transcript = currentanno$ensembl,
                gene = currentanno$symbol, strand = currentstrand,
                window = windowvec, coord = windcoordvec)
            return(res)}, expbed, windcoordvec, nbwindows)
    stopCluster(cl)
    nbwindcheck <- unique(sapply(windflist, nrow))
    if (!isTRUE(all.equal(length(nbwindcheck), 1)) ||
        !isTRUE(all.equal(nbwindcheck, 200)))
        stop("Problem in the nb of windows per transcript retrieved")
    windf <- do.call("rbind", windflist)

    return(windf)
}

makewindowsbedtools <- function(expbed, nbwindows, nbcputrans, verbose = TRUE) {

    ## Filtering out intervals smaller than nbwindows
    idxsmall <- which((expbed$end - expbed$start) < nbwindows)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        message("Excluding ", lsmall, "/", nrow(expbed), " annotations that ",
        "are too short.")
        expbed <- expbed[-idxsmall,]
    }

    ## Splitting each transcript into "nbwindows" windows
    if (verbose) message("\t Splitting ", nrow(expbed), " transcript into ",
        nbwindows, " windows data.frame")
    windcoordvec <- seq_len(nbwindows)
    winddf <- .divideannoinwindows(expbed, windcoordvec, nbwindows,
        nbcputrans)

    return(winddf)
}


!!!!!!!!!!!!!!!!! functions moved to backup folder in file wmeanrelatedfunctionspreprocess.R



##################
# MAIN
##################

createfolder(robjoutputfold)

## Reading the information about experiments
exptab <- read.csv(exptabpath, header = TRUE)
checkexptab(exptab)

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

## Keeping "transcript" annotations
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

## Combine the annotations
message("Combine the annotations")
protcodbed <- cbind(protcodbed, biotype = "protein-coding")
lncrnabed <- cbind(lncrnabed, biotype = "lncRNA")
allannobed <- rbind(protcodbed, lncrnabed)

saveRDS(allannobed, file.path(robjoutputfold, "allannobed.rds"))

## Make windows for all annotations
message("Make windows for all annotations")
idxpar <- grep("PAR_Y", allannobed$ensembl)
if (!isTRUE(all.equal(length(idxpar), 0)))
    allannobed <- allannobed[-idxpar, ]

allwindowsbed <- makewindowsbedtools(expbed = allannobed, nbwindows = windsize,
    nbcputrans = nbcputrans)
saveRDS(allwindowsbed, file.path(robjoutputfold, "allwindowsbed.rds"))


## Retrieving the values of the bedgraph files, removing black lists and keeping
## high mappability scores
message("Removing the black list and keeping scores with high mappability")
if (is.null(blacklistshpath)) {
    message("Retrieving the black list online")
    blacklistgr <- createblacklist(blacklistname, outputfolder)
} else {
    blacklistbed <- read.delim(blacklistshpath, header = FALSE)
}

maptrackbed <- read.delim(maptrackpath, header = FALSE)
saveRDS(maptrackbed, file.path(robjoutputfold, "maptrackbed.rds"))

expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
!!bedgraphwmeanlist <- retrieveandfilterfrombg(exptab, blacklistgr,
    maptrackgr, nbcpubg, expnamevec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Creating a rowid that will be used for merging
message("Adding rowid for each bedgraph")
bedgraphwmeanlist <- mclapply(bedgraphwmeanlist, function(tab) {
    rowidvec <- paste(tab$biotype, tab$seqnames, tab$start, tab$end, tab$strand,
        tab$gene, tab$transcript, paste0("frame", tab$window),
        paste0("coord", tab$coord), sep = "_")
    ## Inserting rowid col after transcript
    res <- data.frame(seqnames = tab$seqnames, start = tab$start, end = tab$end,
        strand = tab$strand, gene = tab$gene, biotype = tab$biotype,
        window = tab$window, coord = tab$coord, transcript = tab$transcript,
        rowid = rowidvec, tab[, grep("score", colnames(tab))])
    colnames(res)[ncol(res)] <- colnames(tab)[grep("score", colnames(tab))]
    return(res)
}, mc.cores = nbcpubg)

message("Joining the elements of each bedgraph")
start_time <- Sys.time()
completeframedf <- purrr::reduce(bedgraphwmeanlist, dplyr::full_join,
    by = c("seqnames", "start", "end", "strand", "gene", "biotype", "window",
    "coord", "transcript", "rowid"))
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time)
saveRDS(completeframedf, file = "/g/romebioinfo/tmp/preprocessing/completeframedf-blacklistfile.rds") # nolint

saveRDS(completeframedf, file = "/g/romebioinfo/tmp/preprocessing/completeframedf.rds") # nolint
