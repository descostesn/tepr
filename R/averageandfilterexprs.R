#' Calculate Average Expression and Filter Transcript Data
#'
#' This function calculates the average expression levels for transcripts from
#' a provided expression data frame and filters out transcripts based on a
#' specified expression threshold. The function also renames the columns in the
#' output data frame to include mean expression values.
#'
#' @param expdf A data frame containing expression data that should have
#'              columns named 'condition', 'replicate', and 'strand'.
#' @param alldf A data frame containing all transcript-related information,
#'              including biotype, chromosome, coordinates, transcript, gene,
#'              strand, window, ID and scores retrieved from the bedgraph
#'              files.
#' @param expthres A numeric value specifying the expression threshold.
#'                 Transcripts with average expression values below this
#'                 threshold will be filtered out from the returned transcript
#'                 vector.
#' @param verbose A logical value indicating whether to print progress messages
#'                 (default is FALSE).
#'
#' @return A list containing:
#'         \item{maintable}{The original data frame containing all transcript
#'          data.}
#'         \item{exptranstab}{A character vector of transcripts that meet the
#'                            filtering criteria.}
#'
#' @examples
#' # Example usage of averageandfilterexprs
#' result <- averageandfilterexprs(expdf, alldf, expthres = 10, verbose = TRUE)
#'
#' @importFrom dplyr group_by summarize filter select bind_rows arrange pull
#' @importFrom tidyselect all_of contains
#'
#' @export

averageandfilterexprs <- function(expdf, alldf, expthres, verbose = FALSE) { # nolint

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
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint

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
                dplyr::filter(strand == strandname) %>% # nolint
                dplyr::select(gene, transcript, strand, # nolint
                tidyselect::contains(directname))  %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ !is.na(.))) %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ . > expthres))
            return(dfstrand)
        }, unique(dfbytranscript$strand), unique(expdf$strand),
            MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)

    exptranstab <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>%
        dplyr::arrange(transcript) %>% dplyr::pull(transcript) # nolint

    return(list(maintable = alldf, exptranlist = exptranstab))
}
