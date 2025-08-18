.retrievebgval <- function(currentpath, currentchrom, chromlength, showmemory,
    verbose) {

    if (verbose) message("\t\t\t Reading ", currentpath)
    whichchrom <- GenomicRanges::GRanges(paste0(currentchrom, ":1-",
        chromlength))
    valgr <- rtracklayer::import.bedGraph(currentpath, which = whichchrom)
    if (verbose) message("\t\t\t Converting to tibble")
    valdf <- as.data.frame(valgr)
    colnames(valdf) <- c("chrom", "start", "end", "width", "strand", "score")
    valtib <- tibble::as_tibble(valdf)
    rm(valdf)
    if (showmemory) print(gc()) else invisible(gc())
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
    res <- list(currenttrans, idxscore)
    return(res)
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
            stop("\n\t The number of frames should be equal to ",
                "windsize: ", windsize, " for transcript ",
                unique(currenttrans$transcript.window), ".\n")
        idxscorereplace <- match(dupframenbvec,
            currenttrans$window.window)

        if (!isTRUE(all.equal(dupframenbvec,
            currenttrans$window.window[idxscorereplace])))
                stop("\n\t Problem in replacing scores by wmean, ",
                    "contact the developer.\n")
            currenttrans[idxscorereplace, "score"] <- wmeanvec
    }
    return(currenttrans)
}

.retrievemaptrack <- function(maptrackpath, showtime, showmemory, currentchrom,
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
            suppressWarnings(maptrackbedchrom <- rtracklayer::import(
                maptrackbedfile, which = whichchrom))
            maptrackbedchrom <- as.data.frame(maptrackbedchrom)

            idxvec <- match(c("name", "width"), colnames(maptrackbedchrom))
            if (!isTRUE(all.equal(length(idxvec), 2)))
                stop("\n\t Columns 'name' and 'width' were not found in the ",
                    "maptrack file. This should not happen. Contact the ",
                    "developer.\n")
            maptrackbedchrom <- maptrackbedchrom[, -idxvec]
            colnames(maptrackbedchrom) <- c("chrom", "start", "end", "id",
                "mapscore")
            maptracktib <- tibble::as_tibble(maptrackbedchrom)

            rm(filename, maptrackbedfile, whichchrom, maptrackbedchrom)
            if (showmemory) print(gc()) else invisible(gc())

            if (!is.na(saveobjectpath)) {
                if (verbose) message("\t\t Saving mappability track to ",
                    maptrackbedobjfile)
                saveRDS(maptracktib, maptrackbedobjfile)
            }
        } else {
            if (verbose) message("\t\t Loading mappability track from ",
                    "existing rds object")
            maptracktib <- readRDS(maptrackbedobjfile)
        }

        if (showtime) {
            end_time_maptrackreading <- Sys.time()
            timing <- end_time_maptrackreading - start_time_maptrackreading
            message("\t\t ## Read maptrack: ", format(timing, digits = 2))
        }

        return(maptracktib)
}

.retrievechromlength <- function(chromtab, currentchrom) {

    idxchrom <- which(GenomeInfoDb::seqnames(chromtab) == currentchrom)
    chromlength <- as.numeric(GenomeInfoDb::seqlengths(chromtab)[idxchrom])
    return(chromlength)
}

#' Retrieve chromosome lengths and information for a specified genome.
#'
#' @description
#' This function connects to the UCSC Genome Browser database using the
#' `rtracklayer` package to retrieve chromosome information. It returns a
#' `Seqinfo` object, filtering out unwanted chromosomes such as mitochondrial
#' DNA (`chrM`) and those with alternative contigs (indicated by an underscore
#' `_`).
#'
#' @usage
#' retrievechrom(genomename, verbose, filterchrom = TRUE)
#' 
#' @param genomename A character string specifying the UCSC genome name (e.g.,
#' "hg19" or "mm10").
#' @param verbose A logical value. If `TRUE`, the function will print messages
#'   during execution, including a list of the chromosomes being kept.
#' @param filterchrom A logical value. If `TRUE`, mitochondrial and non-canonical
#'  chromosomes are removed. Default is \code{TRUE}.
#'
#' @return A `Seqinfo` object containing the names and lengths of the main
#'   chromosomes for the specified genome.
#'
#' @examples
#' # This example requires an internet connection to the UCSC database
#' hg19_chroms <- retrievechrom(genomename = "hg19", verbose = TRUE)
#' hg19_chroms
#'
#' @importFrom rtracklayer SeqinfoForUCSCGenome
#' @importFrom GenomeInfoDb seqnames
#' @export
retrievechrom <- function(genomename, verbose, filterchrom = TRUE) {

    if (verbose) message("Retrieving chromosome lengths")
    chromtab <- rtracklayer::SeqinfoForUCSCGenome(genomename)
    if (is.null(chromtab))
        stop("\n\t The genome ", genomename, " was not found with the function",
        " rtracklayer::SeqinfoForUCSCGenome. Check the spelling or verify",
        " if the genome is available on UCSC. The connection to UCSC can ",
        "also have some hickup. You can callagain the function using the ",
        "chromtab parameter: chromtab <- rtracklayer::SeqinfoForUCSCGenome(",
        "genomename).\n")
    if (filterchrom) {
        idxkeep <- GenomeInfoDb::seqnames(chromtab)[grep("_|chrM",
        GenomeInfoDb::seqnames(chromtab), perl = TRUE, invert = TRUE)]
        chromtab <- chromtab[idxkeep,]
    }
    if (verbose) message("\t Working on: ",
        paste(GenomeInfoDb::seqnames(chromtab), collapse="/"))
    return(chromtab)
}

.retrieveannoscores <- function(currentstrand, allwindchromtib, valtib, # nolint
    showmemory, verbose) {

        ## Declaration to tackle CMD check
        strand <- NULL

        ## Keeping information on the correct strand
        if (verbose) message("\t\t Retrieving information on strand ", # nolint
            currentstrand)
        if (isTRUE(all.equal(currentstrand, "plus")))
            retrievedstrand <- "+"
        else
            retrievedstrand <- "-"
        allwindstrand <- allwindchromtib %>%
            dplyr::filter(strand == as.character(retrievedstrand)) # nolint

        ## Retrieving scores on annotations of strand
        if (verbose) message("\t\t Retrieving scores on annotations of strand") # nolint
        suppressWarnings(annoscores <- valr::bed_intersect(valtib,
            allwindstrand, suffix = c("", ".window")))

        rm(valtib, allwindchromtib, allwindstrand)
        if (showmemory) print(gc()) else invisible(gc())
        return(annoscores)
}

.rowidandcols <- function(bytranslist, currentcond, currentrep, # nolint
    currentdirection, showmemory, verbose) {

        ## Declaration to tackle CMD check
        biotype <- chrom <- NULL

        ## Combining transcripts in one table
        if (verbose) message("\t\t Combining transcripts in one table")
        res <- do.call("rbind", bytranslist)
        rm(bytranslist)
        if (showmemory) print(gc()) else invisible(gc())

        if (verbose) message("\t\t Formatting and adding rowid column")
        ## Create rowid string
        rowidvec <- paste(res$transcript, res$gene, res$strand,
            res$window, sep = "_")
        ## Inserting rowid col after window
        res <- res %>% tibble::add_column(rowid = rowidvec,
            .after = "window")
        ## Move biotype col before chrom col
        res <- res %>% dplyr::relocate(biotype, .before = chrom) # nolint
        ## Retrieving score column position
        idxcolscore <- grep("_score", colnames(res))
        ## Creating experiment columns
        expcol <- paste0(currentcond, "_rep", currentrep, ".",
            currentdirection)
        expcolvec <- rep(expcol, nrow(res))
        tmpres <- cbind(res[, -idxcolscore], expcolvec)
        res <- cbind(tmpres, res[, idxcolscore])
        res <- tibble::as_tibble(res)

        rm(tmpres)
        if (showmemory) print(gc()) else invisible(gc())

        return(res)
}
