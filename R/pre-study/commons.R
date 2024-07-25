bedtogr <- function(currentbed, strand = TRUE, symbol = TRUE) {

    grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
            ranges = IRanges::IRanges(start = currentbed[, 2],
                                  end = currentbed[, 3],
                                  names = currentbed[, 4]),
            strand = if (strand) currentbed[, 6] else "+",
            symbol = if (symbol) currentbed[, 5] else NA)
    return(grres)
}

excludeorkeepgrlist <- function(expgr, removegr, removefrom = TRUE,
    ignorestrand = TRUE) {
    ## command retrieved with HelloRanges:
    # nolint - bedtools_intersect("-a protcod.bed -b hg38-blacklist.v2.bed -v")
    resgr <- IRanges::subsetByOverlaps(expgr, removegr,
        invert = removefrom, ignore.strand = ignorestrand)
    return(resgr)
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
    expnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction)
    scorelist <- parallel::mcmapply(function(bwpath, expname, grintervals,
        verbose) {
        if (verbose) message("\t Retrieving bw values for ", expname)
        meanvec <- .retrievemeanfrombw(grintervals, bwpath, verbose)
        return(meanvec)
    }, expdf$path, expnamevec, MoreArgs = list(grintervals, verbose),
        SIMPLIFY = FALSE, mc.cores = nbcpu)

    ## Building matrix of mean scores
    if (verbose) message("\t Building matrix")
    scoremat <- do.call("cbind", scorelist)
    colnames(scoremat) <- expnamevec
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
    windowvec <- .returnwindowvec(dfintervalsrownames)

    ## Final data.frame
    if (verbose) message("\t Creating the final data.frame")
    df <- data.frame(biotype = grname, chr = dfintervals$seqnames,
    start = dfintervals$start, end = dfintervals$end, transcript = transvec,
    gene = dfintervals$symbol, strand = dfintervals$strand, window = windowvec,
    id = paste(transvec, dfintervals$symbol, dfintervals$strand, windowvec,
        sep = "_"))
    df <- cbind(df, scoremat)
    rownames(df) <- NULL

    return(df)
}
