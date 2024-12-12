.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed) {

    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window")

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

.wmeanvec <- function(dupframenbvec, currenttrans) {

    wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans) {

        ## Selecting all rows having a window equal to nbdup
        allframedf <- currenttrans[which(currenttrans$window.window == nbdup), ]

        ## Retrieving start and end of the window to calculate nb of nt
        windowstart <- unique(allframedf$start.window)
        windowend <- unique(allframedf$end.window)
        if (!isTRUE(all.equal(length(windowstart), 1)) ||
            !isTRUE(all.equal(length(windowend), 1)))
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

    return(wmeanvec)
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


.removeblackandlowmap <- function(currenttrans, blacklisttib, idxscore,
    maptracktib) {

        ## Set scores overlapping black list to NA
        resblack <- valr::bed_intersect(currenttrans, blacklisttib)
        if (!isTRUE(all.equal(nrow(resblack), 0))) {
            strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                currenttrans$end, sep = "-")
            strblack <- paste(resblack$chrom, resblack$start.x, resblack$end.x,
                sep = "-")
            idxblack <- as.vector(na.omit(unique(match(strtransvec, strblack))))
            if (isTRUE(all.equal(length(idxblack), 0)))
                stop("Problem in setting scores overlapping black list to NA.",
                    " This should not happen. Contact the developer.")
            currenttrans[idxblack, idxscore] <- NA
        }

        ## Set scores NOT overlapping high map to NA (i.e. scores overlapping
        ## low mappability intervals)
        resmap <- valr::bed_intersect(currenttrans, maptracktib)
        if (!isTRUE(all.equal(nrow(resmap), 0))) {
            ## Compute strtransvec only if no overlap with black list was found
            if (isTRUE(all.equal(nrow(resblack), 0)))
                strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                    currenttrans$end, sep = "-")
            strmap <- paste(resmap$chrom, resmap$start.x, resmap$end.x,
                sep = "-")
            idxmap <- match(strtransvec, strmap)
            idxtosetna <- which(is.na(idxmap)) ## NOT overlapping high map
            currenttrans[idxtosetna, idxscore] <- NA
        }
        # rm(resblack, resmap, strtransvec, strblack, strmap)
        # invisible(gc())
        return(currenttrans)
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

.meanblackhighbytrans <- function(bgscorebytrans, windsize, currentname,
    blacklisttib, maptracktib, saveobjectpath, nbcputrans, reload, verbose) {

        currentobj <- file.path(saveobjectpath, paste0(currentname,
            "-translist.rds"))
        if (!reload || !file.exists(currentobj)) {
            bytranslist <- parallel::mclapply(bgscorebytrans,
                function(currenttrans, windsize, currentname, blacklisttib,
                    maptracktib) {

                        ## Reordering rows according to window number
                        idxwind <- order(currenttrans$window.window)
                        currenttrans <- currenttrans[idxwind, ]

                        ## Identifying duplicated windows that will be used to
                        ## compute a weighted mean.
                        currenttrans <- .dupidx(currenttrans, windsize)

                        ## Remove columns corresponding to bedgraph, move score
                        ## column at the end, remove the .window suffix from
                        ## column names, Add exp name prefix to column score
                        res <- .formatcurrenttranscols(currenttrans,
                            currentname)
                        currenttrans <- res[[1]]
                        idxscore <- res[[2]]

                        ## Set scores overlapping black list and low map to NA
                        idxchrom <- which(maptracktib$chrom == unique(
                            currenttrans$chrom))
                        maptracktibchrom <- maptracktib[idxchrom, ]
                        currenttrans <- .removeblackandlowmap(currenttrans,
                            blacklisttib, idxscore, maptracktibchrom)

                        return(currenttrans)

            }, windsize, currentname, blacklisttib, maptracktib,
                mc.cores = nbcputrans)
            if (!is.na(saveobjectpath)) {
                if (verbose) message("\t\t\t Saving ", currentobj)
                saveRDS(bytranslist, file = currentobj)
            }
        } else {
            if (verbose) message("\t\t\t Loading ", currentobj)
            bytranslist <- readRDS(currentobj)
        }

        return(bytranslist)
}

.retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed, # nolint
    nbcputrans, allwindowsbed, expnamevec, windsize, saveobjectpath, showtime,
    reload, verbose) {

        if (verbose) message("\t Converting annotations' windows, blacklist,",
            " and mappability track to tibble")
        tibres <- .convertotibble(allwindowsbed, blacklistbed, maptrackbed)
        allwindtib <- tibres[[1]]
        blacklisttib <- tibres[[2]]
        maptracktib <- tibres[[3]]
        if (verbose) message("\t Deleting objects and free memory")
        rm(tibres, allwindowsbed, blacklistbed, maptrackbed)
        invisible(gc())

        ## Looping on each experiment bg file
        if (verbose) message("\t For each bedgraph file")
        bedgraphlistwmean <- mapply(function(currentpath, currentname,
            currentstrand, allwindtib, blacklisttib, maptracktib, windsize,
            nbcputrans, saveobjectpath, verbose, showtime, reload) {

            ## Retrieving bedgraph values
            if (verbose) message("\n\t\t Retrieving begraph values for ",
                currentname)
            invisible(gc())
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

            ## Retrieving scores on annotations of strand
            if (verbose) message("\t\t Retrieving scores on annotations of ",
                "strand")
            suppressWarnings(annoscores <- valr::bed_intersect(valtib,
                allwindstrand, suffix = c("", ".window")))

            ## Splitting the scores by transcript
            if (verbose) message("\t\t Splitting the scores by transcript")
            trsfact <- factor(annoscores$transcript.window)
            bgscorebytrans <- split(annoscores, trsfact)
            if (verbose) message("\t\t Deleting objects and free memory")
            rm(trsfact, valtib)
            invisible(gc())

             ## For each transcript compute the weighted means for each window.
             ## The weight is calculated if a window contains more than one
             ## score
             if (verbose) message("\t\t For each transcript compute the ",
                "weighted means and set scores overlapping black list and low ",
                "mappability to NA. It takes a while.")
             if (showtime) start_time_bytranslist <- Sys.time()
             bytranslist <- .meanblackhighbytrans(bgscorebytrans, windsize,
                currentname, blacklisttib, maptracktib, saveobjectpath,
                nbcputrans, reload, verbose)
            if (showtime) {
                end_time_bytranslist <- Sys.time()
                timing <- end_time_bytranslist - start_time_bytranslist
                message("\t\t ## Exp treated in: ", timing) # nolint
            }

            if (!isTRUE(all.equal(unique(sapply(bytranslist, nrow)), windsize)))
                stop("All elements of the list should contain ", windsize,
                    " rows. This should not happen. Contact the developer.")

            ## Combining transcripts in one table
            if (verbose) message("\t\t Combining transcripts in one table")
            res <- do.call("rbind", bytranslist)

            if (verbose) message("\t\t Deleting objects and free memory")
            rm(bgscorebytrans, bytranslist)
            invisible(gc())
            return(res)

        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, nbcputrans, saveobjectpath,
        verbose, showtime, reload), SIMPLIFY = FALSE)

        if (verbose) message("\t\t Free memory")
        invisible(gc())
        return(bedgraphlistwmean)
}

.retrievemaptrackbed <- function(maptrackpath, showtime, saveobjectpath, reload,
    verbose) {

        if (showtime) start_time_maptrackreading <- Sys.time()
        maptrackbedobjfile <- file.path(saveobjectpath, "maptrackbed.rds")

        if (!reload || !file.exists(maptrackbedobjfile)) {
            if (verbose) message("Reading the mappability track (the file is ",
                "big, be patient)")
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

        if (showtime) {
            end_time_maptrackreading <- Sys.time()
            timing <- end_time_maptrackreading - start_time_maptrackreading
            message("\t\t ## read maptrack in: ", timing) # nolint
        }

        return(maptrackbed)
}

.loadbgprocessing <- function(exptab, blacklistbed, maptrackbed, allwindowsbed,
        windsize, nbcputrans, showtime, saveobjectpath, reload, verbose) {

            if (verbose) message("Removing scores within black list intervals,",
            " keeping those on high mappability regions, and computing ",
            "weighted means.")
            expnamevec <- paste0(exptab$condition, exptab$replicate,
                exptab$direction)
            if (showtime) start_time_bglistwmean <- Sys.time()
            bedgraphlistwmean <- .retrieveandfilterfrombg(exptab, blacklistbed,
                maptrackbed, nbcputrans, allwindowsbed, expnamevec, windsize,
                saveobjectpath, showtime, reload, verbose)
            if (showtime) {
                end_time_bglistwmean <- Sys.time()
                timing <- end_time_bglistwmean - start_time_bglistwmean
                message("\t\t ## Built bedgraphlistwmean in: ", timing) # nolint
            }

            if (!is.na(saveobjectpath)) {
                if (verbose) message("Saving bedgraphlistwmean as an rds ",
                    "object")
                saveRDS(bedgraphlistwmean, file.path(saveobjectpath,
                    "bedgraphlistwmean.rds"))
            }

            return(bedgraphlistwmean)
}


#' Process Bedgraph Files by Removing Blacklist Scores and Keeping High
#' Mappability Scores
#'
#' @description
#' This function processes bedgraph files by filtering out scores overlapping
#' blacklisted regions and retaining those in high-mappability regions. It
#' computes weighted means for scores within overlapping windows and supports
#' parallel processing for efficiency.
#'
#' @usage
#' blacklisthighmap(maptrackpath, blacklistshpath, exptabpath, nbcputrans,
#'  allwindowsbed, windsize, saveobjectpath = NA, reload = FALSE,
#'  showtime = FALSE, verbose = TRUE)
#'
#' @param maptrackpath Character string. Path to the mappability track file.
#' @param blacklistshpath Character string. Path to the blacklist regions file.
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param nbcputrans Number of CPU cores to use for transcript-level operations.
#' @param allwindowsbed Data frame. BED-formatted data frame obtained with the
#'  function makewindows.
#' @param windsize Window size for splitting transcripts into intervals.
#' @param saveobjectpath Path to save intermediate R objects. Default is `NA`
#'  and R objects are not saved.
#' @param reload Logical. If `TRUE`, reloads existing saved objects to avoid
#'  recomputation. Default is `FALSE`. If the function failed during object
#'  saving, make sure to delete the corresponding object.
#' @param showtime Logical. Whether to display timing information.
#' @param verbose Logical. Whether to display detailed progress messages.
#'
#' @return A list of data frames where each entry corresponds to the processed
#'  scores for an experiment. Scores outside high-mappability regions or in
#'  blacklisted regions are set to `NA`.
#'
#' @details
#' The function involves the following steps:
#' 1. Reading and converting the blacklist, mappability track, and annotation
#'  windows into tibbles.
#' 2. Reading experiment metadata to identify bedgraph files.
#' 3. For each bedgraph file:
#'    - Extracting scores based on strand information.
#'    - Filtering scores overlapping blacklisted regions or outside
#'  high-mappability intervals.
#'    - Computing weighted means for overlapping windows.
#' 4. Combining processed scores into a single data frame for each experiment.
#'
#' The function can save intermediate objects to disk and reload them for
#' efficiency in repeated runs.
#'
#' @examples
#' # Define paths to required files
#' maptrackpath <- "path/to/maptrack.bed"
#' blacklistshpath <- "path/to/blacklist.bed"
#' exptabpath <- "path/to/experiments.csv"
#' allwindowsbed <- data.frame(...)
#'
#' # Run the function
#' results <- blacklisthighmap(
#'     maptrackpath = maptrackpath,
#'     blacklistshpath = blacklistshpath,
#'     exptabpath = exptabpath,
#'     nbcputrans = 4,
#'     allwindowsbed = allwindowsbed,
#'     windsize = 100,
#'     saveobjectpath = "output/",
#'     reload = FALSE,
#'     showtime = TRUE,
#'     verbose = TRUE
#' )
#'
#' @seealso
#' [makewindows]
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr relocate filter
#' @importFrom rtracklayer import.bedGraph
#' @importFrom valr bed_intersect
#'
#' @export

blacklisthighmap <- function(maptrackpath, blacklistshpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, saveobjectpath = NA, reload = FALSE,
    showtime = FALSE, verbose = TRUE) {

        if (showtime) start_time_fun <- Sys.time()

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- read.csv(exptabpath, header = TRUE)

        if (verbose) message("Reading the black list")
        blacklistbed <- read.delim(blacklistshpath, header = FALSE)

        ## For the mappability track, reading can be skept by loading the object
        ## if it exists
        maptrackbed <- .retrievemaptrackbed(maptrackpath, showtime,
            saveobjectpath, reload, verbose)

        ## Removing scores within black list intervals, keeping those on high
        ## mappability regions, and computing weighted means.
        bedgraphlistwmean <- .loadbgprocessing(exptab, blacklistbed,
            maptrackbed, allwindowsbed, windsize, nbcputrans, showtime,
            saveobjectpath, reload, verbose)

        if (showtime) {
            end_time_fun <- Sys.time()
            timing <- end_time_fun - start_time_fun
            message("\t\t ## All bedgraphs processed in: ", timing) # nolint
        }

        return(bedgraphlistwmean)
}
