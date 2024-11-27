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
            !!.windsizevec <- function()
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
