.retrievekneeandmax <- function(condvec, transtable) { # nolint

reslist <- lapply(condvec, function(cond, transtable) {
    difffxname <- paste0("diff_Fx_", cond)
    difffxvec <- transtable[, difffxname]
    ## If equality of difference within the same gene it takes the closest
    ## knee from the TSS # nolint
    resrow <- transtable[which(difffxvec == max(difffxvec)), ] %>% # nolint
        dplyr::slice_min(.data$coord, n = 1) # nolint
    res <- data.frame(resrow$coord, resrow[, difffxname])
    colnames(res) <- c(paste0("knee_AUC_", cond), paste0("max_", difffxname))
    return(res)
}, transtable)

return(reslist)
}

#' Identify the Knee and Max ECDF Differences for Each Transcript
#'
#' @description
#' This function identifies the knee point (i.e., point of maximum change) and
#' the maximum difference in the empirical cumulative distribution function
#' (ECDF) for each transcript, across different experimental conditions.
#'
#' @usage
#' kneeid(transdflist, expdf, nbcpu = 1, showtime = FALSE, verbose = TRUE)
#'
#' @param transdflist A list of data frames where each data frame contains
#'    transcript data with ECDF values for each condition.
 #' @param expdf A data frame containing experimental information including a
 #'   \code{condition} column.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'  parallel computation. The parallelization is performed on the elements of
#'  transdflist. Defaults to 1.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return A data frame where each row corresponds to a transcript and contains
#'  the coordinates of the knee point and the maximum ECDF difference for each
#'  condition.
#'
#' @examples
#' # Assuming transdflist is a list of transcript data frames and expdf contains
#' # conditions for each experiment:
#' # result <- kneeid(transdflist, expdf, nbcpu = 4, verbose = TRUE)
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr slice_min
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export

kneeid <- function(transdflist, expdf, nbcpu = 1, showtime = FALSE,
  verbose = TRUE) {

  if (showtime) start_time <- Sys.time()
  condvec <- unique(expdf$condition)
  bytransres <- parallel::mclapply(transdflist, function(transtable, condvec) {
      bycondreslist <- .retrievekneeandmax(condvec, transtable)
       return(cbind(transcript = transtable$transcript[1],
        do.call("cbind", bycondreslist)))
    }, condvec, mc.cores = nbcpu)
  res <- do.call("rbind", bytransres)

  if (showtime) {
      end_time <- Sys.time()
      message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
  }

  return(res)
}
