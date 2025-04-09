#' Count NA values per transcript and condition
#'
#' @description
#' This function takes a list of expression data frames, a condition
#' information data frame, and counts the number of NA values for each
#' transcript based on strand and condition. NA represent missing scores that
#' were filtered out from the black list and mappability track. The function
#' operates in parallel on transcripts to speed up the process using multiple
#' CPU cores.
#'
#' @usage
#' countna(allexprsdfs, expdf, nbcpu = 1, showtime = FALSE, verbose = TRUE)
#'
#' @param allexprsdfs A list of data frames containing expression data. The
#'  first element is assumed to be the main table. The second element is a
#'  vector of transcript names that passed the filtering of
#'  'averageandfilterexprs'.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'  parallel computation on transcripts. The number of transcripts is equal to
#'  the number of lines provided as input of 'averageandfilterexprs'.
#'  Defaults to \code{1}.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return A data frame where each row corresponds to a transcript, along with
#'  its associated gene, strand, and the count of NA values.
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Calculating averageandfilterexprs to call countNA
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
#'         showtime = FALSE, verbose = FALSE)
#'
#' ## Testing countna
#' res <- countna(avfilttest, expdf, nbcpu = 1, verbose = FALSE)
#'
#'
#' @importFrom parallel mclapply
#' @importFrom stats na.omit
#'
#' @seealso
#' [averageandfilterexprs]
#'
#' @export

countna <- function(allexprsdfs, expdf, nbcpu = 1, showtime = FALSE,
  verbose = TRUE) {

  if (showtime) start_time <- Sys.time()
  if (verbose) message("\n\t ## Counting NA values") # nolint
  maintable <- allexprsdfs[[1]]
  scorecolvec <- grep("_score", colnames(maintable), value = TRUE)
  condvec <- unique(expdf$condition)

  ## Splitting the table by each transcript
  if (verbose) message("\t Splitting the table by each transcript") # nolint
  transdflist <- split(maintable, factor(maintable$transcript))

  ## For each transcript
  nabytranslist <- parallel::mclapply(transdflist,
    function(transtable, scorecolvec, condvec) {
        ## Filters the score columns according to the strand of the transcript
        str <- .extractstr(transtable)
        colnamestr <- scorecolvec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]
        ## If there is only one condition and one replicate
        if (isTRUE(all.equal(length(colnamestr), 1))) {
          scoremat <- as.matrix(scoremat)
          colnames(scoremat) <- colnamestr
        }

        ## Counting NA for each condition (c=condition, m=matrix, n=colnames)
        res <- sapply(condvec, function(c, m, n) {
          idxcol <- grep(c, n)
          if (!isTRUE(all.equal(length(idxcol), 1))) {
            return(length(which(apply(m[, idxcol], 1,
              function(x) all(is.na(x))))))
          } else { ## Only one replicate
            return(length(which(is.na(m[, idxcol]))))
          }
        }, scoremat, colnamestr)

        ## Retrieving total NA and transcript info
        if (!isTRUE(all.equal(length(condvec), 1)) && 
          !isTRUE(all.equal(res[[1]], res[[2]])))
            stop("\n\t Number of NA is different between conditions for ",
              unique(transtable$gene), ": ", res[[1]], " - ", res[[2]],
              ". This should not happen. Contact the developer.\n")
        ## I drop the other NA columns because it is the same value for all the
        ## conditions (NA depends on blacklist and unmmapable region)
        countna <- res[1]
        info <- data.frame(gene = unique(transtable$gene),
          transcript = unique(transtable$transcript),
          strand = unique(transtable$strand))
        resmat <- cbind(info, Count_NA = countna)

        return(resmat)
    }, scorecolvec, condvec, mc.cores = nbcpu)

  finalres <- do.call("rbind", nabytranslist)

  if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
  }

  return(finalres)
}
