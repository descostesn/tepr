.computewindflist <- function(expbed, nbwindows) {

    n <- nrow(expbed)

    ## Compute window sizes for all genes at once (vectorized)
    lgene <- expbed$end - expbed$start
    windowsize <- round(lgene / nbwindows)
    missingbp <- lgene %% nbwindows

    ## Create index vectors: each gene gets nbwindows rows
    gene_idx <- rep(seq_len(n), each = nbwindows)
    window_num <- rep(seq_len(nbwindows), times = n)

    ## Vectorized coordinate computation
    ws <- windowsize[gene_idx]
    st <- expbed$start[gene_idx]

    coor1 <- st + (window_num - 1L) * ws
    coor2 <- st + window_num * ws

    ## Add remainder bp to the last window of each gene
    is_last <- window_num == nbwindows
    coor2[is_last] <- coor2[is_last] + missingbp[gene_idx[is_last]]

    ## Build the result data.frame in one shot
    windf <- data.frame(biotype = expbed$biotype[gene_idx],
        chr = expbed$chrom[gene_idx], coor1 = coor1,
        coor2 = coor2, transcript = expbed$ensembl[gene_idx],
        gene = expbed$symbol[gene_idx], strand = expbed$strand[gene_idx],
        window = window_num, stringsAsFactors = FALSE)

    return(windf)
}


.divideannoinwindows <- function(expbed, nbwindows) {

    windf <- .computewindflist(expbed, nbwindows)

    ## Validate: total rows should be nrow(expbed) * nbwindows
    expected_rows <- nrow(expbed) * nbwindows
    if (!isTRUE(all.equal(nrow(windf), expected_rows)))
        stop("\n[tepr] Error: Incorrect window count per transcript.\n",
            "  Contact the developer.\n")

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
    winddf <- .divideannoinwindows(expbed, nbwindows)

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
#'  specified number of windows (`windsize`). It uses vectorized operations to
#' efficiently split transcripts into fixed-size windows. The result includes
#' metadata for each window, such as its chromosome, start and end coordinates,
#' associated gene, and the window number.
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
#' @seealso
#' [retrieveanno]
#'
#' @export

makewindows <- function(allannobed, windsize, nbcputrans = 1, verbose = TRUE,
    saveobjectpath = NA, showtime = FALSE) {

        if (showtime) start_time <- Sys.time()
        
        ## Making windows for all annotations
        if (verbose) message("Making windows for all annotations")

        # Filtering out annotations with "PAR_Y" in the ensembl name
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
