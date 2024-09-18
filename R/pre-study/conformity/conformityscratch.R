library("rtracklayer")
library("GenomicRanges")
library("tibble")
library("valr")


##################
# PARAMETERS
##################

bgvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score"

allbgnicpath <- "/g/romebioinfo/tmp/preprocessing/backup/bedgraphwmeanlist.rds"
allwindowspath <- "/g/romebioinfo/tmp/preprocessing/allwindowsbed.rds"

exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv" # nolint
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint

nbcpubg <- 1
windsize <- 200

##################
#FUNCTIONS
##################



##################
# MAIN
##################

## Reading all windows bed
allwindowsbed <- readRDS(allwindowspath)

## Reading exptab, black list, and maptrack
exptab <- read.csv(exptabpath, header = TRUE)
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
blacklistbed <- read.delim(blacklistshpath, header = FALSE)
maptrackbed <- read.delim(maptrackpath, header = FALSE)


##################### TEST
## This is the ctrl rep1 fwd
bgvic <- read.delim(bgvicpath, header = FALSE)

## Selecting ctrl rep1 fwd
allbgnic <- readRDS(allbgnicpath)
names(allbgnic) <- gsub(".bg","",basename(names(allbgnic)))
bgnic <- allbgnic[["ctrl_rep1.forward"]]
rm(allbgnic)
gc()
## Selecting the lines corresponding to the gene ARF5
bgvicarf <- bgvic[which(bgvic$V6 == "ARF5"), ]
bgnicarf <- bgvic[which(bgnic$gene == "ARF5"), ]
allwindarf <- allwindowsbed[which(allwindowsbed$gene == "ARF5"), ]
#allwindowsgr <- bedtogr(allwindarf, allwindows = TRUE)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed, verbose) {

    if (verbose) message("Converting annotations' windows to tibble") # nolint
    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    if (verbose) message("Converting blacklist to tibble") # nolint
    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

    if (verbose) message("Converting mappability track to tibble") # nolint
    colnames(maptrackbed) <- c("chrom", "start", "end", "id", "mapscore")
    maptracktib <- tibble::as_tibble(maptrackbed)

    return(list(allwindtib, blacklisttib, maptracktib))
}

.retrievebgval <- function(currentpath, verbose) {

    valgr <- rtracklayer::import.bedGraph(currentpath)
    if (verbose) message("\t\t Converting to tibble")
    valdf <- as.data.frame(valgr)
    colnames(valdf) <- c("chrom", "start", "end", "width", "strand", "score")
    valtib <- tibble::as_tibble(valdf)
    return(valtib)
}

.removeblacklist <- function(allwindstrand, valtib, currentstrand,
    blacklisttib, verbose) {

        if (verbose) message("\t Retrieving scores on annotations of strand ", # nolint
                currentstrand)
        suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
                suffix = c("", ".window")))

        ## Removing black list
        if (verbose) message("\t Keeping scores outside blacklist intervals") # nolint
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)

        return(resblack)
}

.retrieveonhighmap <- function(resblack, maptracktib, currentchrom) {

    ## Keeping scores on high mappability track
    resmap <-  valr::bed_intersect(resblack,
        maptracktib  %>% dplyr::filter(chrom == currentchrom), # nolint
        suffix = c(".bg", ".maphigh"))
    colnames(resmap) <- gsub(".window.bg", ".window", colnames(resmap))

    ## Removing mapping columns and duplicates
    resmap <- resmap[, -grep(".maphigh|.overlap|.source", colnames(resmap))]
    resmap <- resmap %>% dplyr::distinct(chrom, start.bg, end.bg, # nolint
        start.window, end.window, .keep_all = TRUE) # nolint
    invisible(gc())
    return(resmap)
}

.uniqueformatcolnames <- function(currenttrans) {

    ## Verifying uniformity of chrom, transcript, and genes
    uniquechrom <- as.character(unique(currenttrans$chrom))
    uniquetrans <- as.character(unique(currenttrans$transcript.window))
    uniquegene <- as.character(unique(currenttrans$gene.window))

    if (!isTRUE(all.equal(length(uniquechrom), 1)) ||
        !isTRUE(all.equal(length(uniquetrans), 1)) ||
        !isTRUE(all.equal(length(uniquegene), 1)))
            stop("chrom, transcript, and genes should be unique, this should", # nolint
                " not happen. Contact the developper.") # nolint

    ## Renaming window and coord columns removing the suffix
    colnamevec <- colnames(currenttrans)
    colnames(currenttrans)[which(colnamevec == "window.window")] <- "window"
    colnames(currenttrans)[which(colnamevec == "coord.window")] <- "coord"
    colnames(currenttrans)[which(colnamevec == "gene.window")] <- "gene"
    colnames(currenttrans)[which(colnamevec == "transcript.window")] <- "transcript" # nolint

    return(list(currenttrans, uniquechrom, uniquetrans, uniquegene))

}

.retrievemissingwind <- function(idxnavec, allwindstrand, currenttrans,
    uniquechrom, uniquetrans, uniquegene) {

    ## For each missing window whose number is contained in idxnavec
    missingrowslist <- lapply(idxnavec, function(idxna, allwindstrand) {

        ## Retrieving the line of the missing window in allwindstrand
        idxmissing <-  which(allwindstrand$chrom == uniquechrom &
            allwindstrand$transcript == uniquetrans &
            allwindstrand$gene == uniquegene &
            allwindstrand$window == idxna)

        if (!isTRUE(all.equal(length(idxmissing), 1)))
            stop("Problem in retrieving the missing window, this should not ",
                "happen. Contact the developper.")

        ## Below the bedgraph information columns are set to NA. These columns will be removed later # nolint
        ## The score is set to NA since it is a missing value resulting from removing black list and low mappability (keeping high mappability) # nolint
        ## Filling the other columns with the line retrieved in allwindstrand # nolint
        windstrandrow <- allwindstrand[idxmissing, ]
        resmissing <- data.frame(chrom = windstrandrow$chrom,
            start.bg = NA, end.bg = NA, width.bg = NA, strand.bg = "*", ## Set the bedgraph info # nolint
            score.bg = NA, ## Set the score to NA to keep track of missing values # nolint
            biotype.window = windstrandrow$biotype,
            start.window = windstrandrow$start,
            end.window = windstrandrow$end,
            transcript = windstrandrow$transcript, gene = windstrandrow$gene,
            strand.window = windstrandrow$strand, window = windstrandrow$window,
            coord = windstrandrow$coord)

        return(resmissing)
        }, allwindstrand)

    missingrowsdf <- do.call("rbind", missingrowslist)
    currenttrans <- rbind(currenttrans, missingrowsdf)
    currenttrans <- currenttrans[order(currenttrans$coord), ]

    return(currenttrans)
}

.arrangewindows <- function(currenttrans, windsize, allwindstrand) {

    res <- .uniqueformatcolnames(currenttrans)
    currenttrans <- res[[1]]
    uniquechrom <- res[[2]]
    uniquetrans <- res[[3]]
    uniquegene <- res[[4]]

    ## Identifying the missing window in currenttrans
    idx <- match(seq_len(windsize), unique(currenttrans$window))
    idxnavec <- which(is.na(idx))

    ## If some windows are missing
    if (!isTRUE(all.equal(length(idxnavec), 0)))
        currenttrans <- .retrievemissingwind(idxnavec,
            allwindstrand, currenttrans, uniquechrom,
            uniquetrans, uniquegene)

    scorebg <- currenttrans$score.bg
    idxremovebg <- grep(".bg", colnames(currenttrans))
    currenttrans <- currenttrans[, -idxremovebg]
    currenttrans <- cbind(currenttrans, scorebg)
    idxscore <- which(colnames(currenttrans) == "scorebg")
    scorename <- paste0(currentname, "_score") # nolint
    colnames(currenttrans)[idxscore] <- scorename
 
    return(currenttrans)
}

retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed,
    nbcpubg, allwindowsbed, expnamevec, windsize, verbose = TRUE) {

        tibres <- .convertotibble(allwindowsbed, blacklistbed, maptrackbed,
            verbose)
        allwindtib <- tibres[[1]]
        blacklisttib <- tibres[[2]]
        maptracktib <- tibres[[3]]

        ## Looping on each experiment bg file
        bedgraphlist <- parallel::mcmapply(function(currentpath, currentname,
            currentstrand, allwindtib, blacklisttib, maptracktib, nbcpuchrom,
            windsize, verbose) {

            ## Retrieving bedgraph values
            if (verbose) message("\t Retrieving values for ", currentname) # nolint
            valtib <- .retrievebgval(currentpath, verbose)

            ## Keeping window coordinates on the correct strand
            allwindstrand <- allwindtib %>%
                dplyr::filter(strand == as.character(currentstrand)) # nolint

            ## Overlapping scores with anno on correct strand and remove
            ## blacklist
            resblack <- .removeblacklist(allwindstrand, valtib, currentstrand,
                blacklisttib, verbose)

            ## Processing by chromosomes because of size limits, the mappability
            ## track has too many rows
            chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
            resmaplist <- lapply(chromvec, function(currentchrom, allwindstrand,
                currentname, resblack, maptracktib) {

                if (verbose) message("\t\t over ", currentchrom)

                if (verbose) message("\t\t\t Keeping scores on high ",
                    "mappability track")
                resmap <- .retrieveonhighmap(resblack, maptracktib,
                    currentchrom)

                ## Processing data per transcript
                message("\t\t\t Building scoring results by transcript")
                bgscorebytrans <- split(resmap,
                    factor(resmap$transcript.window))

                message("\t\t\t Setting missing windows scores to NA")
                #currenttrans=bgscorebytrans[[1]]
                bytranslist <- lapply(bgscorebytrans,
                    function(currenttrans, windsize, allwindstrand,
                        currentname) {

                            currenttrans <- .arrangewindows(currenttrans,
                                windsize, allwindstrand)

                }, windsize, allwindstrand, currentname)

            }, allwindstrand, currentname, resblack, maptracktib)
        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, verbose), SIMPLIFY = FALSE,
        mc.cores = nbcpubg)
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



