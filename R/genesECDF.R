.coordandfilter <- function(str, transtable, nbrows) { # nolint

  if (isTRUE(all.equal(str, "minus"))) {
    coordvec <- seq(from = nbrows, to = 1, by = -1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable[order(coordvec), ]
    transtable <- transtable %>% dplyr::select(!dplyr::matches("plus"))
  } else {
    coordvec <- seq(from = 1, to = nbrows, by = 1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable %>% dplyr::select(!dplyr::matches("minus"))
  }
  return(transtable)
}

.computeecdf <- function(transtable, rounding, nbrows) { # nolint

        ## Declaration to tackle CMD check
        variable <- NULL

        ## Retrieving keyword plus or minus
        str <- .extractstr(transtable)

        ## Create the coordinate column and select scores having the righ
        ## orientation
        transtable <- .coordandfilter(str, transtable, nbrows)
        colnamevec <- colnames(transtable)
        colscorevec <- colnamevec[grep("_score", colnamevec)]

        ## Filling the NA of the score columns of the right strand with
        ## tidyr::fill in the downup direction
        transtable <- transtable %>% tidyr::fill(tidyr::contains("score"),
           .direction = "downup")

        ## Computing ecdf
        suppressWarnings(dflong <- transtable %>%
            tidyr::gather(key = "variable", value = "value", colscorevec))
        dflong[, "value_round"] <- round(dflong$value * rounding)
        ecdflist <- lapply(unique(dflong$variable), function(currentvar) {
            dfsubset <- subset(dflong,
              subset = variable == currentvar)
            dfexpanded <- dfsubset[rep(seq_len(nrow(dfsubset)),
                dfsubset$value_round), ]
            funecdf <- stats::ecdf(dfexpanded[, "coord"])
            dfsubset$Fx <- funecdf(dfsubset$coord)
            return(dfsubset)
        })
        resecdf <- dplyr::bind_rows(ecdflist)

        ## Shrink the results back to the transtable keeping ecdf columns
        res <- dplyr::select(tidyr::pivot_wider(resecdf,
          names_from = "variable", values_from = c("value", "value_round",
          "Fx")), -tidyselect::contains("value_round"))

        ## Removing strand from column names
        colnames(res) <- gsub(paste0(".", str), "", colnames(res))

        return(res)
}

#' Compute ECDF for Genes Based on Expression Data
#'
#' @description
#' This function calculates the empirical cumulative distribution function
#' (ECDF) for expressed genes across multiple transcripts. It processes the
#' expression data to filter out non-expressed transcripts, compute ECDF values
#' for each transcript, and combine the results into a unified data frame. The
#' function operates in parallel for speed optimization.
#'
#' @usage
#' genesECDF(allexprsdfs, nbcpu = 1, rounding = 10,
#' showtime = FALSE, verbose = TRUE)
#'
#' @param allexprsdfs A list of data frames where the first element is the main
#'    expression data frame and the second element contains the names of the
#'    expressed transcripts (see 'averageandfilterexprs').
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'    parallel computation. Default is \code{1}.
#' @param rounding An integer specifying the rounding factor for computing ECDF.
#'    Default is \code{10}.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'    Default is \code{TRUE}.
#'
#' @return A list containing two elements:
#' \item{concatdf}{A data frame with ECDF results for each transcript.}
#' \item{nbrows}{An integer indicating the number of rows in each transcript
#'  table.}
#'
#' @details
#' The function performs several steps:
#' \enumerate{
#'   \item Filters the main expression table to retain only the expressed
#'      transcripts.
#'   \item Splits the data by each transcript.
#'   \item For each transcript, computes ECDF values for the score columns
#'      while respecting the strand orientation ("plus" or "minus").
#'   \item Combines the ECDF results into a final data frame.
#' }
#'
#' The function uses parallel processing to compute ECDF for each transcript
#'  simultaneously, making it faster on systems with multiple CPU cores.
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Calculating averageandfilterexprs and countNA to call genesECDF
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
#'         showtime = FALSE, verbose = FALSE)
#' countnatest <- countna(avfilttest, expdf, nbcpu = 1, verbose = FALSE)
#'
#' ## Testing genesECDF 
#' resecdf <- genesECDF(avfilttest, verbose = FALSE)
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @importFrom tidyr fill pivot_wider
#' @importFrom tidyselect contains
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stats ecdf
#'
#' @seealso
#' [averageandfilterexprs]
#'
#' @export

genesECDF <- function(allexprsdfs, nbcpu = 1, rounding = 10, # nolint
  showtime = FALSE, verbose = TRUE) {

    if (showtime) start_time <- Sys.time()
    if (verbose) message("\n\t ## Computing ecdf")
    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]
    maincolnamevec <- colnames(maintable)

    ## Filtering the main table to keep only the expressed transcripts
    if (verbose) message("\t Filtering to keep only the expressed transcripts") # nolint
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    if (isTRUE(all.equal(length(idxnoexpr), 0))) {
      if (verbose)
        warning("[tepr] Warning: All transcripts are expressed.",
          immediate. = TRUE)
    } else {
      maintable <- maintable[-idxnoexpr, ]
    }

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    if (verbose) message("\t Splitting the table by each transcript") # nolint
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
    ecdflist <- parallel::mclapply(transdflist, function(transtable,
        rounding, nbrows, maincolnamevec) {

        res <- .computeecdf(transtable, rounding, nbrows)
        return(res)
    }, rounding, nbrows, maincolnamevec, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }

    finalres <- list(concatdf, nbrows)
    return(finalres)
}
