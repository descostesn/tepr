countna <- function(allexprsdfs, expdf, nbcpu, verbose = FALSE) {

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
        if (!isTRUE(all.equal(res[[1]], res[[2]])))
            stop("Number of NA is different between conditions. This should ",
                "not happen. Contact the developer.")
        ## I drop the other NA columns because it is the same value for all the
        ## conditions (NA depends on blacklist and unmmapable region)
        countna <- res[1]
        info <- data.frame(gene = unique(transtable$gene),
          transcript = unique(transtable$transcript),
          strand = unique(transtable$strand))
        return(cbind(info, Count_NA = countna))
    }, scorecolvec, condvec, mc.cores = nbcpu)

  return(do.call("rbind", nabytranslist))
}
