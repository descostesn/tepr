.condcolidx <- function(currentcond, df) {
    idxcond <- grep(currentcond, colnames(df))
    if (isTRUE(all.equal(length(idxcond), 0)))
        stop("\n[tepr] Error: Condition not found.\n",
            "  Condition '", currentcond, "' missing in column names.\n",
            "  Ensure same expdf used in all functions. Contact developer.\n")
    return(idxcond)
}

.idxscorefx <- function(df, idxcond) {
    idxcondfx <- grep("Fx", colnames(df[idxcond]))
    idxcondval <- grep("value_", colnames(df[idxcond]))
    if (isTRUE(all.equal(length(idxcondfx), 0)) ||
        isTRUE(all.equal(length(idxcondval), 0)))
        stop("\n[tepr] Error: Missing columns in meandifference.\n",
            "  'Fx' or 'value_' columns not found. Contact developer.\n")
    idxcondlist <- list(value = idxcond[idxcondval],
            Fx = idxcond[idxcondfx])
    return(idxcondlist)
}

.creatematdiff <- function(condvec, resmean) {

  categoryvec <- c("value", "Fx")
  matdifflist <- lapply(categoryvec, function(currentcat, condvec, resmean) {
    meancolnames <- paste("mean", currentcat, condvec, sep = "_")

    ## Generating all combinations of elements (combn not good)
    idxvec <- seq_len(length(condvec))
    matidx <- matrix(c(idxvec, rev(idxvec)), ncol = 2)

    ## Generating differences of columns
    difflist <- apply(matidx, 2, function(idxvec, meancolnames, resmean,
        currentcat, condvec) {
          ## The function rowDiffs of the package matrixStats substracts the second argument to the first one. To respect the code just above, # nolint
          ## The indexes must be inverted with rev: meancolnames[rev(idxvec)]] -> for instance, given the two columns "mean_value_HS" and "mean_value_ctrl" # nolint
          ## as input, the function rowDiffs will do the subtraction "mean_value_ctrl" - "mean_value_HS" # nolint
          res <- matrixStats::rowDiffs(as.matrix(
            resmean[, meancolnames[rev(idxvec)]]))
          colnamestr <- paste("Diff", paste0("mean", currentcat),
            paste(condvec[idxvec], collapse = "_"), sep = "_")
          res <- as.vector(res)
          attr(res, "name") <- colnamestr
          return(res)
      }, meancolnames, resmean, currentcat, condvec, simplify = FALSE)

      ## Combining vectors into a matrix and defining col names
      diffmat <- do.call("cbind", difflist)
      colnames(diffmat) <- sapply(difflist, function(x) attributes(x)$name)
      return(diffmat)
  }, condvec, resmean)

  ## Building a matrix from the diff on values and Fx
  matdiff <- do.call("cbind", matdifflist)
  return(matdiff)
}

.meandiffscorefx <- function(idxcondlist, df, nbwindows, currentcond,
    colnamevec, verbose) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, nbwindows,
            currentcond, colnamevec, verbose) {
            if (verbose) {
              message("\t Calculating average and difference between ",
                "replicates for columns '", idxname, "' of ", currentcond)
              if (isTRUE(all.equal(length(idxvalvec), 1)))
                warning("[tepr] Warning: Only one replicate. ",
                  "Copying scores to mean columns.", immediate. = TRUE)
            }

            ## Calculating the column of mean scores for currentcond
            ## The result is a data.frame made of a single column
            if (length(idxvalvec) >= 2) {
                meandf <- data.frame(rowMeans(df[, idxvalvec], na.rm = FALSE))
            } else {
                meandf <- as.data.frame(df[, idxvalvec])
            }
            colnames(meandf) <- paste0("mean_", idxname, "_", currentcond)

            if (isTRUE(all.equal(idxname, "Fx"))) {
                diffres <- meandf - df$coord / nbwindows
                colnames(diffres) <- paste0("diff_", idxname, "_", currentcond)
                res <- cbind(meandf, diffres)
            } else {
                res <- meandf
            }
            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, nbwindows,
            currentcond, colnamevec, verbose), SIMPLIFY = FALSE)

        return(meandifflist)
}

#' Compute Mean and Differences of Scores for Each Condition
#'
#' @description
#' This function calculates the mean values, mean Fx (ECDF) and ECDF differences
#' (Fx) for expression data, across different experimental conditions. If only
#' one condition is provided, skips computation of mean differences.
#'
#' @usage
#' meandifference(resultsecdf, expdf, nbwindows, showtime = FALSE,
#' verbose = TRUE)
#'
#' @param resultsecdf A data frame containing ECDF results for each transcript
#'  and condition (see genesECDF).
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param nbwindows An integer representing the number of windows (or segments)
#'  in each transcript.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return A data frame that contains, for each condition:
#' \itemize{
#'   \item Mean values for the "value" and "Fx" columns (e.g.,
#'    \code{mean_value_ctrl}, \code{mean_Fx_ctrl}).
#'   \item Differences between the \code{Fx} column and coordinate ratios
#'    (e.g., \code{diff_Fx_ctrl}).
#' }
#' If only one condition is provided, the differences on mean columns are not
#' performed.
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Calculating necessary results
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' avfilt <- averageandfilterexprs(expdf, transdf, expthres,
#'         showtime = FALSE, verbose = FALSE)
#' countna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#'
#' ## Testing meandifference
#' meandifftest <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#'
#' @importFrom dplyr bind_rows
#' @importFrom matrixStats rowDiffs
#'
#' @export

meandifference <- function(resultsecdf, expdf, nbwindows, showtime = FALSE,
  verbose = TRUE) {

    if (showtime) start_time <- Sys.time()
    if (verbose) message("\n\t ## Computing meandifference")
    ## for each condition, creates three columns:
    ##   - "mean_value_ctrl", "mean_Fx_ctrl", "diff_Fx_ctrl"
    ##   - "mean_value_HS", "mean_Fx_HS", "diff_Fx_HS"
    condvec <- unique(expdf$condition)
    rescondlist <- lapply(condvec, function(currentcond, df, nbwindows,
      verbose) {

        if (verbose) message("\t Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)
        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxscorefx(df, idxcond)

        ## The difference is used to calculate the AUC later on
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, nbwindows,
          currentcond, colnamevec, verbose)
        names(meandifflist) <- NULL

        meandiffres <- do.call("cbind", meandifflist)
        return(meandiffres)
    }, resultsecdf, nbwindows, verbose)
    resmean <- do.call("cbind", rescondlist)

    ## Computing all differences on mean columns
    if (!isTRUE(all.equal(length(condvec), 1))) {

      if (verbose) message("\t Computing all differences on mean columns")
      matdiff <- .creatematdiff(condvec, resmean)

      res <- cbind(resmean, matdiff)
      if (!isTRUE(all.equal(nrow(resultsecdf), nrow(res))))
          stop("\n[tepr] Error: Row count mismatch.\n",
              "  Mean/diff results differ from ecdf results.\n",
              "  If expdf has 2 conditions, contact developer.\n",
              "  Otherwise use teprmulti(). Check with showallcomp(expdf).\n")
    } else {
      if (verbose) message("\t There is only one condition. Skip Computing all",
        " differences on mean columns.")
      res <- resmean
    }

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }

    res <- cbind(resultsecdf, res)
    return(res)
}
