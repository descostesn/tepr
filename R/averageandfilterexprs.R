#' Calculate Average Expression and Filter Transcript Data
#'
#' @description
#' This function calculates the average expression levels for transcripts from
#' a provided expression data frame and filters out transcripts based on a
#' specified expression threshold. The function also renames the columns in the
#' output data frame to include mean expression values.
#'
#' @usage
#' averageandfilterexprs(expdf, alldf, expthres, showtime = FALSE,
#' verbose = TRUE)
#'
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param alldf A data frame containing all transcript-related information,
#'              including biotype, chromosome, coordinates, transcript, gene,
#'              strand, window, ID and scores retrieved from the bedgraph
#'              files.
#' @param expthres A numeric value specifying the expression threshold.
#'                 Transcripts with average expression values below this
#'                 threshold will be filtered out from the returned transcript
#'                 vector.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical value indicating whether to print progress messages
#'                 Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'         \item{maintable}{The original data frame containing all transcript
#'          data.}
#'         \item{exptranstab}{A character vector of transcripts that meet the
#'                            filtering criteria.}
#'
#' @examples
#' # Example usage of averageandfilterexprs
#' # result <- averageandfilterexprs(expdf, alldf, expthres = 10)
#'
#' @importFrom dplyr group_by summarize filter select bind_rows arrange pull
#' @importFrom tidyselect all_of contains
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export

averageandfilterexprs <- function(expdf, alldf, expthres, showtime = FALSE, # nolint
    verbose = TRUE) {

    if (showtime) start_time <- Sys.time()
    ## Verify the conformity of the experiment table
    checkexptab(expdf)

    ## Adding column names to alldf
    infocolnames <- c("biotype", "chr", "coor1", "coor2", "transcript",
        "gene", "strand", "window", "id")
    expcolnames <- unlist(apply(expdf, 1, function(x) {
        res <- paste0(x["condition"], "_rep", x["replicate"], ".", x["strand"])
        return(c(res, paste(res, "score", sep = "_")))
    }, simplify = FALSE))
    colnames(alldf) <- c(infocolnames, expcolnames)

    scorecolvec <- expcolnames[grep("_score", expcolnames)]

    ## Calculate the average expression per transcript (over each frame)
    if (verbose) message("\t Calculating average expression per transcript") # nolint
    dfbytranscript <- alldf %>% dplyr::group_by(.data$transcript) %>%
        dplyr::summarize(gene = .data$gene[1],
            strand = .data$strand[1],
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean"))

    ## Remove a line if it contains only values < expthres (separating strands)
    if (verbose) message("\t Removing lines with values < expthres") # nolint
    dfstrandlist <- mapply(function(strandname, directname, dfbytrans,
        expthres) {
            if ((isTRUE(all.equal(strandname, "plus")) &&
                isTRUE(all.equal(directname, "reverse"))) ||
                (isTRUE(all.equal(strandname, "minus")) &&
                isTRUE(all.equal(directname, "forward"))))
                    stop("Strand and direction do not match, contact the ",
                        "developer")
            dfstrand <- dfbytranscript %>%
                dplyr::filter(.data$strand == strandname) %>%
                dplyr::select(gene, transcript, strand,
                tidyselect::contains(directname))  %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ !is.na(.))) %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ . > expthres))
            return(dfstrand)
        }, unique(dfbytranscript$strand), unique(expdf$strand),
            MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)

    exptranstab <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>%
        dplyr::arrange(.data$transcript) %>%
        dplyr::pull(.data$transcript)

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t ## Analysis performed in: ", format(timing, digits = 2))
    }
    return(list(maintable = alldf, exptranlist = exptranstab))
}
