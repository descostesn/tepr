bedtogr <- function(currentbed, strand = TRUE) {

    grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
            ranges = IRanges::IRanges(start = currentbed[, 2],
                                  end = currentbed[, 3],
                                  names = currentbed[, 4]),
            strand = if (strand) currentbed[, 6] else "+")
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
    ## command retrieved with HelloRanges:
    ## bedtools_makewindows("-n 200 -b stdin.bed") # nolint
    ## Note: In R, bedtools does not have the "-i srcwinnum" option
    res <- GenomicRanges::tile(expgr, n = binsize)
    res <- unlist(res, use.names = FALSE)

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
