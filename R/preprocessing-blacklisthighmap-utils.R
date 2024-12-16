.retrievebgval <- function(currentpath, verbose) {

    if (verbose) message("\t\t Reading ", currentpath)
    valgr <- rtracklayer::import.bedGraph(currentpath)
    if (verbose) message("\t\t Converting to tibble")
    valdf <- as.data.frame(valgr)
    colnames(valdf) <- c("chrom", "start", "end", "width", "strand", "score")
    valtib <- tibble::as_tibble(valdf)
    return(valtib)
}

.formatcurrenttranscols <- function(currenttrans, currentname) {

    ## Remove columns corresponding to bedgraph
    idxcolbg <- match(c("start", "end", "width", "strand", ".overlap"),
        colnames(currenttrans))
    currenttrans <- currenttrans[, -idxcolbg]

    ## Move score column at the end
    currenttrans <- dplyr::relocate(currenttrans, "score",
        .after = "window.window")

    ## Remove the .window suffix from column names
    colnames(currenttrans) <- gsub(".window", "", colnames(currenttrans))

    ## Add exp name prefix to column score
    idxscore <- which(colnames(currenttrans) == "score")
    colnames(currenttrans)[idxscore] <- paste(currentname, "score", sep = "_")

    return(list(currenttrans, idxscore))
}

.dupidx <- function(currenttrans, windsize) {

    dupidx <- which(duplicated(currenttrans$window.window))
    if (!isTRUE(all.equal(length(dupidx), 0))) {

        ## Computing the weighted mean for each duplicated window
        dupframenbvec <- unique(currenttrans$window.window[dupidx])
        wmeanvec <- .wmeanvec(dupframenbvec, currenttrans)

        ## Remove duplicated frames and replace scores by wmean
        currenttrans <- currenttrans[-dupidx, ]
        if (!isTRUE(all.equal(nrow(currenttrans), windsize)))
            stop("The number of frames should be equal to ",
                "windsize: ", windsize, " for transcript ",
                unique(currenttrans$transcript.window))
        idxscorereplace <- match(dupframenbvec,
            currenttrans$window.window)

        if (!isTRUE(all.equal(dupframenbvec,
            currenttrans$window.window[idxscorereplace])))
                stop("Problem in replacing scores by wmean, ",
                    "contact the developer.")
            currenttrans[idxscorereplace, "score"] <- wmeanvec
    }
    return(currenttrans)
}

.retrievemaptrackbed <- function(maptrackpath, showtime, currentchrom,
    chromlength, saveobjectpath, reload, verbose) {

        if (showtime) start_time_maptrackreading <- Sys.time()
        filename <- paste0("maptrackbed-", currentchrom, ".rds")
        maptrackbedobjfile <- file.path(saveobjectpath, filename)

        if (!reload || !file.exists(maptrackbedobjfile)) {
            if (verbose) message("\t\t Reading the mappability track on the",
                " specified chromosome")
            ## Reading file on chrom and converting to data.frame
            maptrackbedfile <- rtracklayer::BEDFile(maptrackpath)
            whichchrom <- GenomicRanges::GRanges(
                paste0(currentchrom, ":1-", chromlength))
            maptrackbedchrom <- rtracklayer::import(maptrackbedfile,
                which = whichchrom)
            maptrackbedchrom <- as.data.frame(maptrackbedchrom)

            idxvec <- match(c("name", "width"), colnames(maptrackbedchrom))
            if (!isTRUE(all.equal(length(idxvec), 2)))
                stop("Columns 'name' and 'width' were not found in the ",
                    "maptrack file. This should not happen. Contact the ",
                    "developer.")
            maptrackbedchrom <- maptrackbedchrom[, -idxvec]
            colnames(maptrackbedchrom) <- c("chrom", "start", "end", "id",
                "mapscore")
            maptracktib <- tibble::as_tibble(maptrackbedchrom)

            rm(filename, maptrackbedfile, whichchrom, maptrackbedchrom)
            invisible(gc())

            if (!is.na(saveobjectpath)) {
                if (verbose) message("Saving mappability track to ",
                    maptrackbedobjfile)
                saveRDS(maptracktib, maptrackbedobjfile)
            }
        } else {
            if (verbose) message("Loading mappability track from existing rds ",
                    "object")
            maptracktib <- readRDS(maptrackbedobjfile)
        }

        if (showtime) {
            end_time_maptrackreading <- Sys.time()
            timing <- end_time_maptrackreading - start_time_maptrackreading
            message("\t\t ## read maptrack in: ", timing) # nolint
        }

        return(maptracktib)
}

.retrievechromlength <- function(chromtab, currentchrom) {

    idxchrom <- which(GenomeInfoDb::seqnames(chromtab) == currentchrom)
    chromlength <- as.numeric(GenomeInfoDb::seqlengths(chromtab)[idxchrom])
    return(chromlength)
}

.retrievechrom <- function(genomename, verbose) {

    if (verbose) message("Retrieving chromosome lengths")
    chromtab <- rtracklayer::SeqinfoForUCSCGenome(genomename)
    if (is.null(chromtab))
        stop("The genome ", genomename, " was not found with the function ",
        " rtracklayer::SeqinfoForUCSCGenome. Check the spelling or verify",
        " if the genome is available on UCSC.")
    idxkeep <- GenomeInfoDb::seqnames(chromtab)[grep("_|chrM",
        GenomeInfoDb::seqnames(chromtab), perl = TRUE, invert = TRUE)]
    chromtab <- chromtab[idxkeep,]
    if (verbose) message("\t Working on: ",
        paste(GenomeInfoDb::seqnames(chromtab), collapse="/"))
    return(chromtab)
}
