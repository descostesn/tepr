.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname, # nolint
                " should be unique, contact the developer.") # nolint
}

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

.computeecdf <- function(transtable, expdf, rounding, nbrows) { # nolint

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
            dfsubset <- subset(dflong, subset = variable == currentvar) # nolint
            dfexpanded <- dfsubset[rep(seq_len(nrow(dfsubset)),
                dfsubset$value_round), ]
            funecdf <- ecdf(dfexpanded[, "coord"])
            dfsubset$Fx <- funecdf(dfsubset$coord)
            return(dfsubset)
        })
        resecdf <- dplyr::bind_rows(ecdflist)

        ## Shrink the results back to the transtable keeping ecdf columns
        res <- resecdf %>% tidyr::pivot_wider(.,
            names_from = "variable",
            values_from = c("value", "value_round", "Fx")) %>%
            dplyr::select(., -tidyselect::contains("value_round"))

        ## Removing strand from column names
        res <- res %>% dplyr::rename_with(~gsub(paste0(".", str), "", .),
                tidyselect::contains(paste0(".", str)))

        return(res)
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1, # nolint
  verbose = FALSE) {

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]
    maincolnamevec <- colnames(maintable)

    ## Filtering the main table to keep only the expressed transcripts
    if (verbose) message("\t Filtering to keep only the expressed transcripts") # nolint
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    if (isTRUE(all.equal(length(idxnoexpr), 0)))
      warning("All the transcripts are expressed", immediate. = TRUE) # nolint
    else
      maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    if (verbose) message("\t Splitting the table by each transcript") # nolint
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        rounding, nbrows, maincolnamevec) {

        res <- .computeecdf(transtable, expdf, rounding, nbrows)
        return(res)
    }, expdf, rounding, nbrows, maincolnamevec, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(list(concatdf, nbrows))
}
