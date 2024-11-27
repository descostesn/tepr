## For different strings provided in the vector "valvec", perform the filtering
## of gentab on each string using the result of the previous filtering
.grepsequential <- function(valvec, gentab, verbose, invert = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ] ## This line enables sequential grep
    }))
    return(gentab)
}

.sortedbedformat <- function(gencode) {
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


retrieveanno <- function(exptabpath, gencodepath, saveobjectpath = NA,
    verbose = TRUE) {

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Reading the information about experiments
    if (verbose) message("Reading the information about experiments")
    exptab <- read.csv(exptabpath, header = TRUE)
    checkexptab(exptab) # nolint

    ## Reading gencode file
    if (verbose) message("Reading gencode file and filtering")
    gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

    ## Keeping "transcript" annotations
    if (verbose) message("\t Keeping 'transcript' annotations")
    gencode <- gencode[which(gencode$V3 == "transcript"), ]

    ## Selecting Ensembl_canonical transcripts i.e. most representative
    ## transcript of the protein coding gene. This will be the MANE_Select
    ## transcript if there is one, or a transcript chosen by an Ensembl
    ## algorithm otherwise.
    if (verbose) message("\t Selecting Ensembl_canonical transcripts and ",
        "sorting")
    gencodeprotcod <- .grepsequential("MANE_Select", gencode, verbose)
    protcodbed <- .sortedbedformat(gencodeprotcod)

    ## Selecting long non-coding transcripts
    if (verbose) message("\t Selecting long non-coding transcripts and ",
        "sorting")
    lncrna <- .grepsequential(c("lncRNA", "Ensembl_canonical"), gencode,
        verbose)
    removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
    lncrna <- .grepsequential(removevec, lncrna, verbose, invert = TRUE)
    lncrnabed <- .sortedbedformat(lncrna)

    ## Combine the annotations
    if (verbose) message("\t Combine the annotations")
    protcodbed <- cbind(protcodbed, biotype = "protein-coding")
    lncrnabed <- cbind(lncrnabed, biotype = "lncRNA")
    allannobed <- rbind(protcodbed, lncrnabed)

    if (!is.na(saveobjectpath))
        saveRDS(allannobed, file.path(saveobjectpath, "allannobed.rds"))

    return(allannobed)
}


###########################

.windsizevec <- function(currentstart, currentend, nbwindows) {

    lgene <- currentend - currentstart
    windowsize <- round(lgene / nbwindows)
    missingbp <- lgene %% nbwindows
    windsizevec <- rep(windowsize, nbwindows)

    ## Add the missing nb of bp (that is ignore by tile) in the last
    ## element of windsizevec
    if (!isTRUE(all.equal(missingbp, 0)))
        windsizevec[nbwindows] <- windsizevec[nbwindows] + missingbp

    return(windsizevec)
}

.buildveccumsum <- function(currentstart, currentend, nbwindows, currentstrand,
    windcoordvec) {

        ## Compute the vector with the size of each window
        windsizevec <- .windsizevec(currentstart, currentend, nbwindows)

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

        return(list(startvec, endvec, windowvec))
}

.computewindflist <- function(nbcputrans, expbed, windcoordvec, nbwindows) {

    cl <- parallel::makeCluster(nbcputrans)
    windflist <- parallel::parLapply(cl, seq_len(nrow(expbed)),
    function(i, expbed, windcoordvec, nbwindows) {

        ## Retrieve the necessary gene information
        currentanno <- expbed[i, ]
        currentstart <- currentanno$start
        currentend <- currentanno$end
        currentstrand <- currentanno$strand
        windowvec <- windcoordvec

        ## Compute the vector with the size of each window
        ## Building the start and end vectors using the cummulative sum
        ## Inverting start, end, and window vectors if strand is negative
        res <- .buildveccumsum(currentstart, currentend, nbwindows,
            currentstrand, windcoordvec)
        startvec <- res[[1]]
        endvec <- res[[2]]
        windowvec <- res[[3]]

        ## Build the result data.frame containing the coordinates of each
        ## frame alongside window and coord numbers
        res <- data.frame(biotype = currentanno$biotype,
            chr = currentanno$chrom, coor1 = startvec,
            coor2 = endvec,  transcript = currentanno$ensembl,
            gene = currentanno$symbol, strand = currentstrand,
            window = windowvec, coord = windcoordvec)

        return(res)
    }, expbed, windcoordvec, nbwindows)
    parallel::stopCluster(cl)

    return(windflist)
}

.divideannoinwindows <- function(expbed, windcoordvec, nbwindows, nbcputrans) {

    ## Retrieve the necessary gene information
    ## Compute the vector with the size of each window
    ## Building the start and end vectors using the cummulative sum
    ## Inverting start, end, and window vectors if strand is negative
    ## Build the result data.frame containing the coordinates of each
    ## frame alongside window and coord numbers
    windflist <- .computewindflist(nbcputrans, expbed, windcoordvec, nbwindows)

    nbwindcheck <- unique(sapply(windflist, nrow))
    if (!isTRUE(all.equal(length(nbwindcheck), 1)) ||
        !isTRUE(all.equal(nbwindcheck, 200)))
        stop("Problem in the nb of windows per transcript retrieved")
    windf <- do.call("rbind", windflist)

    return(windf)
}

.makewindowsbedtools <- function(expbed, nbwindows, nbcputrans, verbose) {

    ## Filtering out intervals smaller than nbwindows
    idxsmall <- which((expbed$end - expbed$start) < nbwindows)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        message("\t Excluding ", lsmall, "/", nrow(expbed),
            " annotations that are too short.")
        expbed <- expbed[-idxsmall, ]
    }

    ## Splitting each transcript into "nbwindows" windows
    if (verbose) message("\t Splitting ", nrow(expbed), " transcript into ",
        nbwindows, " windows data.frame")
    windcoordvec <- seq_len(nbwindows)
    winddf <- .divideannoinwindows(expbed, windcoordvec, nbwindows,
        nbcputrans)

    return(winddf)
}

makewindows <- function(allannobed, windsize, nbcputrans = 1, verbose = TRUE,
    saveobjectpath = NA) {

    ## Making windows for all annotations
    if (verbose) message("Making windows for all annotations")
    idxpar <- grep("PAR_Y", allannobed$ensembl)
    if (!isTRUE(all.equal(length(idxpar), 0)))
        allannobed <- allannobed[-idxpar, ]
    allwindowsbed <- .makewindowsbedtools(allannobed, windsize, nbcputrans,
        verbose)

    if (!is.na(saveobjectpath))
        saveRDS(allwindowsbed, file.path(saveobjectpath, "allwindowsbed.rds"))

    return(allwindowsbed)
}


######################

.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed) {

    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

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
    blacklisttib) {
        ## Retrieving scores on annotations of strand
        suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
                suffix = c("", ".window")))
        ## Removing black list
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)
        return(resblack)
}

.retrieveonhighmap <- function(resblack, maptracktib, currentchrom) { # nolint

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

.arrangewindows <- function(currenttrans, windsize, allwindstrand,
    currentname) {

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

    idxscore <- which(colnames(currenttrans) == "score.bg")
    scorename <- paste0(currentname, "_score") # nolint
    colnames(currenttrans)[idxscore] <- scorename

    return(list(currenttrans, uniquetrans))
}

.missingandwmean <- function(resmap, windsize, allwindstrand, currentname,
    nbcputrans) {

    ## Splitting the scores kept on the high mappability track by transcript
    bgscorebytrans <- split(resmap, factor(resmap$transcript.window))

    ## Performing identification of missing windows and calculation of weighted
    ## means for each transcript. This is parallelized on nbcputrans CPUs.
    bytranslist <- parallel::mclapply(bgscorebytrans, function(currenttrans,
        windsize, allwindstrand, currentname) {

            ## Identification of missing windows for the current transcript
            ## and setting their scores to NA. Indeed if a window is missing
            ## it is because it was in a black list or a low mappability
            ## interval.
            res <- .arrangewindows(currenttrans, windsize, allwindstrand,
                currentname)
            currenttrans <- res[[1]]
            transname <- res[[2]]

            ## Identifying duplicated windows that will be used to compute
            ## a weighted mean.
            dupidx <- which(duplicated(currenttrans$window))
            colscore <- paste0(currentname, "_score") # nolint

            if (!isTRUE(all.equal(length(dupidx), 0))) {
                dupframenbvec <- unique(currenttrans$window[dupidx])
                ## For each duplicated frame, compute the weighted mean
                wmeanvec <- .computewmeanvec(dupframenbvec, currenttrans,
                    currentname, colscore)
                ## Remove duplicated frames and replace scores by wmean
                currenttrans <- .replaceframeswithwmean(currenttrans, dupidx,
                    windsize, transname, dupframenbvec, colscore, wmeanvec)
             }

             currenttrans <- currenttrans[order(currenttrans$coord), ]
             return(currenttrans)
    }, windsize, allwindstrand, currentname, mc.cores = nbcputrans)

    transdf <- do.call("rbind", bytranslist)
    rm(bytranslist)
    invisible(gc())
    return(transdf)
}

.processingbychrom <- function(maptracktib, allwindstrand, currentname,
    resblack, nbcputrans, windsize, verbose) {

        if (verbose) message("\t\t Retrieving list of chromosomes")
        chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
        if (verbose) message("\t\t Formatting scores")
        bychromlist <- lapply(chromvec, function(currentchrom, allwindstrand,
            currentname, resblack, maptracktib, nbcputrans) {

                if (verbose) message("\t\t\t over ", currentchrom)

                ## Retrieving scores on high mappability intervals
                if (verbose) message("\t\t\t Keeping scores on high ",
                    "mappability track")
                resmap <- .retrieveonhighmap(resblack, maptracktib,
                    currentchrom)

                ## Processing data per transcript for windows and wmean
                message("\t\t\t Setting missing windows scores to NA and",
                    " computing weighted mean for each transcript")
                transdf <- .missingandwmean(resmap, windsize, allwindstrand,
                    currentname, nbcputrans)
                return(transdf)
            }, allwindstrand, currentname, resblack, maptracktib, nbcputrans)

            ## Merging results that were computed on each chromosome
            if (verbose) message("\t\t Merging results that were computed on",
                " each chromome")
            resallchrom <- do.call("rbind", bychromlist)
            rm(bychromlist)
            invisible(gc())
            return(resallchrom)
}

.retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed, # nolint
    nbcputrans, allwindowsbed, expnamevec, windsize, verbose) {

        if (verbose) message("\t Converting annotations' windows, blacklist,",
            " and mappability track to tibble")
        tibres <- .convertotibble(allwindowsbed, blacklistbed, maptrackbed)
        allwindtib <- tibres[[1]]
        blacklisttib <- tibres[[2]]
        maptracktib <- tibres[[3]]

        ## Looping on each experiment bg file
        if (verbose) message("\t For each bedgraph file")
        bedgraphlistwmean <- mapply(function(currentpath, currentname,
            currentstrand, allwindtib, blacklisttib, maptracktib, windsize,
            nbcputrans, verbose) {

            ## Retrieving bedgraph values
            if (verbose) message("\t\t Retrieving begraph values for ",
                currentname)
            valtib <- .retrievebgval(currentpath, verbose)

            ## Keeping window coordinates on the correct strand
            if (verbose) message("\t\t Retrieving coordinates on strand ",
                currentstrand)
            allwindstrand <- allwindtib %>%
                dplyr::filter(strand == as.character(currentstrand)) # nolint

            ## Overlapping scores with anno on correct strand and remove
            ## blacklist
            if (verbose) message("\t\t Keeping scores outside blacklist ",
                "intervals")
            resblack <- .removeblacklist(allwindstrand, valtib, currentstrand,
                blacklisttib)
            rm(valtib)
            invisible(gc())

            ## Processing by chromosomes because of size limits, the mappability
            ## track has too many rows. Formatting scores, keeping those on
            ## high mappability, filling missing windows, and compute wmean
            resallchrom <- .processingbychrom(maptracktib, allwindstrand,
                currentname, resblack, nbcputrans, windsize, verbose)
            return(resallchrom)
        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, nbcputrans, verbose),
        SIMPLIFY = FALSE)

        return(bedgraphlistwmean)
}

## Retrieving the values of the bedgraph files, removing black lists and keeping
## scores landing on high mappability intervals
blacklisthighmap <- function(maptrackpath, blacklistshpath, exptab, nbcputrans,
    allwindowsbed, windsize, saveobjectpath = NA, verbose = TRUE) {

        if (verbose) message("Reading the black list and mappability track")
        blacklistbed <- read.delim(blacklistshpath, header = FALSE)
        maptrackbed <- read.delim(maptrackpath, header = FALSE)

        if (!is.na(saveobjectpath)) {
            if (verbose) message("Saving mappability track as an rds object")
            saveRDS(maptrackbed, file.path(saveobjectpath, "maptrackbed.rds"))
        }

        if (verbose) message("Removing scores within black list intervals, ",
            "keeping those on high mappability regions, and computing ",
            "weighted means.")
        expnamevec <- paste0(exptab$condition, exptab$replicate,
            exptab$direction)
        bedgraphlistwmean <- .retrieveandfilterfrombg(exptab, blacklistbed,
            maptrackbed, nbcputrans, allwindowsbed, expnamevec, windsize,
            verbose)

        return(bedgraphlistwmean)
}
