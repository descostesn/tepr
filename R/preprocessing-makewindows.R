.computewindflist <- function(nbcputrans, expbed, nbwindows) {

    cl <- parallel::makeCluster(nbcputrans)
    windflist <- parallel::parLapply(cl, seq_len(nrow(expbed)),
    function(i, expbed, nbwindows) {

        ## Retrieve the necessary gene information
        currentanno <- expbed[i, ]
        currentstart <- currentanno$start
        currentend <- currentanno$end
        currentstrand <- currentanno$strand

        ## Compute the vector with the size of each window
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
            stop("\n\t Problem in the calculation of windows.\n")

        ## Build the result data.frame containing the coordinates of each
        ## frame alongside window and coord numbers
        res <- data.frame(biotype = currentanno$biotype,
            chr = currentanno$chrom, coor1 = startvec,
            coor2 = endvec,  transcript = currentanno$ensembl,
            gene = currentanno$symbol, strand = currentstrand,
            window = seq_len(nbwindows))

        return(res)
    }, expbed, nbwindows)

    parallel::stopCluster(cl)

    return(windflist)
}


.divideannoinwindows <- function(expbed, nbwindows, nbcputrans) {

    ## Retrieve the necessary gene information
    ## Compute the vector with the size of each window
    ## Building the start and end vectors using the cummulative sum
    ## Inverting start, end, and window vectors if strand is negative
    ## Build the result data.frame containing the coordinates of each
    ## frame alongside window and coord numbers

    windflist <- .computewindflist(nbcputrans, expbed, nbwindows)

    nbwindcheck <- unique(sapply(windflist, nrow))
    if (!isTRUE(all.equal(length(nbwindcheck), 1)) ||
        !isTRUE(all.equal(nbwindcheck, nbwindows)))
        stop("\n\t Problem in the nb of windows per transcript retrieved.",
            " This should not happen. Contact the developer.")
    windf <- do.call("rbind", windflist)

    return(windf)
}

.makewindowsbedtools <- function(expbed, nbwindows, nbcputrans, verbose) {

    ## Filtering out intervals smaller than nbwindows
    idxsmall <- which((expbed$end - expbed$start) < nbwindows)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        if (verbose) message("\t Excluding ", lsmall, "/", nrow(expbed),
            " annotations that are too short.")
        expbed <- expbed[-idxsmall, ]
    }

    ## Splitting each transcript into "nbwindows" windows
    if (verbose) message("\t Splitting ", nrow(expbed), " transcript into ",
        nbwindows, " windows data.frame")
    winddf <- .divideannoinwindows(expbed, nbwindows, nbcputrans)

    return(winddf)
}


#' Split Gene Annotations into Fixed-Size Windows
#'
#' @description
#' This functions uses the annotations filtered from gencode (see retrieveanno).
#' It removes any ensembl names containing "PAR_Y". It filters out intervals
#' smaller than windsize and splits each transcript into "windsize" windows.
#'
#'
#' @usage makewindows(allannobed, windsize, nbcputrans = 1, verbose = TRUE,
#'    saveobjectpath = NA, showtime = FALSE)
#'
#' @param allannobed A data frame which is the result of 'retrieveanno'.
#' @param windsize An integer specifying the number of windows into which each
#'  gene annotation should be divided.
#' @param nbcputrans Number of CPU cores to use for transcript-level operations.
#'  Defaults to 1.
#' @param verbose A logical value indicating whether to display progress
#'  messages. Defaults to `TRUE`.
#' @param saveobjectpath A character string specifying the directory path where
#'  the output object should be saved as an `.rds` file. If `NA`, the object is
#'  not saved. Defaults to `NA`.
#' @param showtime A logical value indicating whether to display the runtime of
#'  the function. Defaults to `FALSE`.
#'
#' @return A data frame containing the split windows for each gene annotation.
#'  The output includes fields such as `biotype`, `chr`, `coor1`, `coor2`,
#'  `transcript`, `gene`, `strand`, and `window`.
#'
#' @details
#' The function filters out annotations with intervals smaller than the
#'  specified number of windows (`windsize`). It uses parallel processing to
#' enhance performance when splitting transcripts into fixed-size windows. The
#' result includes metadata for each window, such as its chromosome, start and
#' end coordinates, associated gene, and the window number.
#'
#' Intermediate functions, such as `.computewindflist` and
#' `.divideannoinwindows`, handle computation and validation of windows. Gene
#' intervals with the "PAR_Y" tag are excluded from the analysis.
#'
#' @examples
#' \donttest{
#' exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
#' gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
#' windsize <- 200
#' 
#' ## Necessary result to call makewindows
#' allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
#'
#' ## Calling makewindows
#' allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)}
#'
#' @importFrom parallel makeCluster parLapply stopCluster
#'
#' @seealso
#' [retrieveanno]
#'
#' @export

makewindows <- function(allannobed, windsize, nbcputrans = 1, verbose = TRUE,
    saveobjectpath = NA, showtime = FALSE) {

        if (showtime) start_time <- Sys.time()
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

        if (showtime) {
            end_time <- Sys.time()
            message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
    }

        return(allwindowsbed)
}
