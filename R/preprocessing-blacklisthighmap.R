.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed) {

    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window")
            # TO REMOVE "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

    colnames(maptrackbed) <- c("chrom", "start", "end", "id", "mapscore")
    maptracktib <- tibble::as_tibble(maptrackbed)

    return(list(allwindtib, blacklisttib, maptracktib))
}

.retrievebgval <- function(currentpath, verbose) {

    if (verbose) message("\t\t Reading ", currentpath)
    valgr <- rtracklayer::import.bedGraph(currentpath)
    if (verbose) message("\t\t Converting to tibble")
    valdf <- as.data.frame(valgr)
    colnames(valdf) <- c("chrom", "start", "end", "width", "strand", "score")
    valtib <- tibble::as_tibble(valdf)
    return(valtib)
}

.removeblacklist <- function(allwindstrand, valtib, blacklisttib) {
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! MOVE OUT THE LINE TO DO MEAN AFTER RETRIEVING BEDGRAPH VALUES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # ## Retrieving scores on annotations of strand
        # suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
        #         suffix = c("", ".window")))
        ## Removing black list
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)
        !!!!!!!!!!! SEE IF NA SHOULD BE ATTRIBUTED HERE
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

    ## Renaming window removing the suffix
    colnamevec <- colnames(currenttrans)
    colnames(currenttrans)[which(colnamevec == "window.window")] <- "window"
    # TO REMOVE colnames(currenttrans)[which(colnamevec == "coord.window")] <- "coord" # nolint
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
            strand.window = windstrandrow$strand, window = windstrandrow$window)
           # TO REMOVE coord = windstrandrow$coord)

        return(resmissing)
        }, allwindstrand)

    missingrowsdf <- do.call("rbind", missingrowslist)
    currenttrans <- rbind(currenttrans, missingrowsdf)
    # TO REMOVE currenttrans <- currenttrans[order(currenttrans$coord), ]

    return(currenttrans)
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

.computewmeanvec <- function(dupframenbvec, currenttrans, currentname,
    colscore) {

    ## For each duplicated frame
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
    wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans, currentname,
        colscore) {

        ## Selecting all rows having a window equal to nbdup
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        allframedf <- currenttrans[which(currenttrans$window == nbdup), ]
!! COMMENTED CODE BECAUSE UNNECESSARY
        # if (isTRUE(all.equal(nrow(allframedf), 1)))
        #     stop("There should be more than one frame selected")

        ## Testing that the start/end of the window is the same for all scores
        ## selected (this should not give an error)
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        windowstart <- unique(allframedf$start.window)
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        windowend <- unique(allframedf$end.window)
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        if (!isTRUE(all.equal(length(windowstart), 1)) ||
            !isTRUE(all.equal(length(windowend), 1)))
                stop("The size of the window is not unique for the frame rows ",
                    "selected, this should not happen, contact the developper.")

        ## Retrieve the nb of overlapping nt for each score
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        overntvec <- apply(allframedf, 1,
            function(x, currentname, windowstart, windowend) {
                nt <- seq(from = x["start.bg"], to = x["end.bg"], by = 1)
                overnt <- length(which(nt >= windowstart & nt <= windowend))
                return(overnt)
            }, currentname, windowstart, windowend)

        ## Computing weighted mean
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        allscores <- as.data.frame(allframedf[,colscore])[[1]]
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        wmean <- weighted.mean(allscores, overntvec)
        return(wmean)
    }, currenttrans, currentname, colscore)
    return(wmeanvec)
}

.replaceframeswithwmean <- function(currenttrans, dupidx, windsize, transname,
    dupframenbvec, colscore, wmeanvec) {

        ## Remove duplicated frames and replace scores by wmean
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        currenttrans <- currenttrans[-dupidx, ]
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        if (!isTRUE(all.equal(nrow(currenttrans), windsize)))
            stop("The number of frames should be equal to windsize: ",
                windsize, " for transcript ", transname)
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        idxscorereplace <- match(dupframenbvec, currenttrans$window)
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        if (!isTRUE(all.equal(dupframenbvec,
            currenttrans$window[idxscorereplace])))
            stop("Problem in replacing scores by wmean, contact the developer.")
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        currenttrans[idxscorereplace, colscore] <- wmeanvec
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
        return(currenttrans)
}

.missingandwmean <- function(resmap, windsize, allwindstrand, currentname,
    nbcputrans) {

!! CODE WAS COPIED AFTER RETRIEVING BG VAL
!!    ## Splitting the scores kept on the high mappability track by transcript
!!    bgscorebytrans <- split(resmap, factor(resmap$transcript.window))

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
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
            dupidx <- which(duplicated(currenttrans$window))
            colscore <- paste0(currentname, "_score") # nolint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! MOVE OUT THE LINE TO DO MEAN AFTER RETRIEVING BEDGRAPH VALUES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! CODE WAS COPIED AFTER RETRIEVING BG VAL
            if (!isTRUE(all.equal(length(dupidx), 0))) {
!! CODE WAS COPIED AFTER RETRIEVING BG VAL
                dupframenbvec <- unique(currenttrans$window[dupidx])
                ## For each duplicated frame, compute the weighted mean
                wmeanvec <- .computewmeanvec(dupframenbvec, currenttrans,
                    currentname, colscore)
                ## Remove duplicated frames and replace scores by wmean
                currenttrans <- .replaceframeswithwmean(currenttrans, dupidx,
                    windsize, transname, dupframenbvec, colscore, wmeanvec)
             }

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             # TO REMOVE currenttrans <- currenttrans[order(currenttrans$coord), ] # nolint
             return(currenttrans)
    }, windsize, allwindstrand, currentname, mc.cores = nbcputrans)

    transdf <- do.call("rbind", bytranslist)
    rm(bytranslist)
    invisible(gc())
    return(transdf)
}

.processingbychrom <- function(maptracktib, allwindstrand, currentname,
    resblack, nbcputrans, windsize, subverbose) {

        if (subverbose) message("\t\t Retrieving list of chromosomes")
        chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
        if (subverbose) message("\t\t Formatting scores")
        bychromlist <- lapply(chromvec, function(currentchrom, allwindstrand,
            currentname, resblack, maptracktib, nbcputrans) {

                if (subverbose) message("\t\t\t over ", currentchrom)

                ## Retrieving scores on high mappability intervals
                if (subverbose) message("\t\t\t Keeping scores on high ",
                    "mappability track")
                resmap <- .retrieveonhighmap(resblack, maptracktib,
                    currentchrom)

                ## Processing data per transcript for windows and wmean
                if (subverbose) message("\t\t\t Setting missing windows scores",
                    " to NA and computing weighted mean for each transcript")
                transdf <- .missingandwmean(resmap, windsize, allwindstrand,
                    currentname, nbcputrans)
                return(transdf)
            }, allwindstrand, currentname, resblack, maptracktib, nbcputrans)

        ## Merging results that were computed on each chromosome
        if (subverbose) message("\t\t Merging results that were computed ",
            "on each chromome")
        resallchrom <- do.call("rbind", bychromlist)
        rm(bychromlist)
        invisible(gc())
        return(resallchrom)
}

.retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed, # nolint
    nbcputrans, allwindowsbed, expnamevec, windsize, verbose, subverbose) {

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

            ## Keeping information on the correct strand
            if (verbose) message("\t\t Retrieving information on strand ",
                currentstrand)
            if (isTRUE(all.equal(currentstrand, "plus")))
                retrievedstrand <- "+"
            else
                retrievedstrand <- "-"
            allwindstrand <- allwindtib %>%
                dplyr::filter(strand == as.character(retrievedstrand)) # nolint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! ASSOCIATE SCORES AND COMPUTE WMEAN HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ## Retrieving scores on annotations of strand
            if (verbose) message("\t\t Retrieving scores on annotations of ",
                "strand")
            suppressWarnings(annoscores <- valr::bed_intersect(valtib,
                allwindstrand, suffix = c("", ".window")))

            ## Splitting the scores by transcript
            if (verbose) message("Splitting the scores by transcript")
            trsfact <- factor(annoscores$transcript.window)
            bgscorebytrans <- split(annoscores, trsfact)

             ## For each transcript compute the weighted means for each window.
             ## The weight is calculated if a window contains more than one
             ## score
             if (verbose) message("For each transcript compute the weighted",
                " mean")
             bytranslist <- parallel::mclapply(bgscorebytrans,
                function(currenttrans, windsize, currentname) {

                    ## Identifying duplicated windows that will be used to
                    ## compute a weighted mean.
                    dupidx <- which(duplicated(currenttrans$window.window))
                    if (!isTRUE(all.equal(length(dupidx), 0))) {
                        dupframenbvec <- unique(currenttrans$window.window[dupidx])

                        ## Computing the weighted mean for each duplicated window
                        wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans) {
                            ## Selecting all rows having a window equal to nbdup
                            allframedf <- currenttrans[which(currenttrans$window.window == nbdup), ]

                            ## Retrieving start and end of the window to calculate nb of nt
                            windowstart <- unique(allframedf$start.window)
                            windowend <- unique(allframedf$end.window)
                            if (!isTRUE(all.equal(length(windowstart), 1)) || !isTRUE(all.equal(length(windowend), 1)))
                                stop("The size of the window is not unique for the frame rows ",
                                    "selected, this should not happen, contact the developper.")

                            ## Retrieve the nb of overlapping nt for each score
                            overntvec <- apply(allframedf, 1, function(x, windowstart, windowend) {
                                nt <- seq(from = x["start"], to = x["end"], by = 1)
                                overnt <- length(which(nt >= windowstart & nt <= windowend))
                                return(overnt)
                            }, windowstart, windowend)

                            ## Computing weighted mean
                            allscores <- as.data.frame(allframedf[, "score"])[[1]]
                            wmean <- weighted.mean(allscores, overntvec)
                            return(wmean)
                        }, currenttrans)

                        ## Remove duplicated frames and replace scores by wmean
                        currenttrans <- currenttrans[-dupidx, ]
                        if (!isTRUE(all.equal(nrow(currenttrans), windsize)))
                            stop("The number of frames should be equal to windsize: ", windsize, " for transcript ", unique(currenttrans$transcript.window))
                        idxscorereplace <- match(dupframenbvec, currenttrans$window.window)

                        if (!isTRUE(all.equal(dupframenbvec, currenttrans$window.window[idxscorereplace])))
                            stop("Problem in replacing scores by wmean, contact the developer.")
                        currenttrans[idxscorereplace, "score"] <- wmeanvec
                    }

                    ## Remove columns corresponding to bedgraph
                    idxcolbg <- match(c("start", "end", "width", "strand", ".overlap"), colnames(currenttrans))
                    currenttrans <- currenttrans[, -idxcolbg]

                    ## Move score column at the end
                    currenttrans <- dplyr::relocate(currenttrans, "score", .after = "window.window")

                    ## Remove the .window suffix from column names
                    colnames(currenttrans) <- gsub(".window", "", colnames(currenttrans))

                    ## Add exp name prefix to column score
                    idxscore <- which(colnames(currenttrans) == "score")
                    colnames(currenttrans)[idxscore] <- paste(currentname, "score", sep = "_")

                        return(currenttrans)
                }, windsize, currentname, mc.cores = nbcputrans)

                if (!isTRUE(all.equal(unique(sapply(bytranslist,nrow)), windsize)))
                    stop("All elements of the list should contain ", windsize, " rows. This should not happen. Contact the developer.")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ## Overlapping scores with anno on correct strand and remove
            ## blacklist
            if (verbose) message("\t\t Keeping scores outside blacklist ",
                "intervals")
            resblack <- .removeblacklist(allwindstrand, valtib, blacklisttib)
            rm(valtib)
            invisible(gc())

            ## Processing by chromosomes because of size limits, the mappability
            ## track has too many rows. Formatting scores, keeping those on
            ## high mappability, filling missing windows, and compute wmean
            if (verbose) message("\t\t Performing operations on each",
                " chromosome. It takes time even with parallelization")
            resallchrom <- .processingbychrom(maptracktib, allwindstrand,
                currentname, resblack, nbcputrans, windsize, subverbose)
            return(resallchrom)
        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, nbcputrans, verbose),
        SIMPLIFY = FALSE)

        return(bedgraphlistwmean)
}

## Retrieving the values of the bedgraph files, removing black lists and keeping
## scores landing on high mappability intervals
blacklisthighmap <- function(maptrackpath, blacklistshpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, saveobjectpath = NA, reload = FALSE,
    verbose = TRUE, subverbose = TRUE) {

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- read.csv(exptabpath, header = TRUE)

        if (verbose) message("Reading the black list")
        blacklistbed <- read.delim(blacklistshpath, header = FALSE)

        ## For the mappability track, reading can be skept by loading the object
        ## if it exists
        maptrackbedobjfile <- file.path(saveobjectpath, "maptrackbed.rds")
        if (!reload || !file.exists(maptrackbedobjfile)) {

            if (verbose) message("Reading the mappability track")
            maptrackbed <- read.delim(maptrackpath, header = FALSE)

            if (!is.na(saveobjectpath)) {
                if (verbose) message("Saving mappability track as an rds ",
                    "object")
                saveRDS(maptrackbed, maptrackbedobjfile)
            }
        } else {
            if (verbose) message("Loading mappability track from existing rds ",
                    "object")
            maptrackbed <- readRDS(maptrackbedobjfile)
        }

        if (verbose) message("Removing scores within black list intervals, ",
            "keeping those on high mappability regions, and computing ",
            "weighted means.")
        expnamevec <- paste0(exptab$condition, exptab$replicate,
            exptab$direction)
        bedgraphlistwmean <- .retrieveandfilterfrombg(exptab, blacklistbed,
            maptrackbed, nbcputrans, allwindowsbed, expnamevec, windsize,
            verbose, subverbose)

        if (!is.na(saveobjectpath)) {
                if (verbose) message("Saving bedgraphlistwmean as an rds ",
                    "object")
                saveRDS(bedgraphlistwmean, file.path(saveobjectpath,
                    "bedgraphlistwmean.rds"))
        }

        return(bedgraphlistwmean)
}
