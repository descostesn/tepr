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
        windsizevec <- .windsizevec(currentstart, currentend, nbwindows)
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

## This functions uses the annotations filtered from gencode (see retrieveanno).
## It removes any ensembl names containing "PAR_Y". It filters out intervals
## smaller than windsize and splits each transcript into "windsize" windows.
makewindows <- function(allannobed, windsize, nbcputrans = 1, verbose = TRUE,
    saveobjectpath = NA) {

    ## Making windows for all annotations
    if (verbose) message("Making windows for all annotations")
    idxpar <- grep("PAR_Y", allannobed$ensembl)
    if (!isTRUE(all.equal(length(idxpar), 0)))
        allannobed <- allannobed[-idxpar, ]
    allwindowsbed <- .makewindowsbedtools(allannobed, windsize, nbcputrans,
        verbose)

    if (!is.na(saveobjectpath)) {
        outfile <- file.path(saveobjectpath, "allwindowsbed.rds")
        if (verbose) message("\t Saving ", outfile)
        saveRDS(allwindowsbed, outfile)
    }

    return(allwindowsbed)
}
