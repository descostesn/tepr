.returninfodf <- function(transtab, nbwindows) { # nolint

    infodf <- transtab  %>%
        dplyr::filter(.data$window == round(nbwindows / 2))  %>%
        dplyr::mutate(window_size = abs(.data$coor2 - .data$coor1), # nolint
          .keep = "all") %>%
        dplyr::select("transcript", "gene", "strand", "window_size") %>%
        dplyr::distinct()

    return(infodf)
}

.checkempty <- function(idx, namestr) {

    if (isTRUE(all.equal(length(idx), 0)))
        stop("\n\t Your condition ", namestr, " was not found in the ",
            "experiment table expdf. Change the parameter controlcondname or",
            " stresscondname.\n")
}

.dauc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1,
    controlcondname = "ctrl", stresscondname = "HS") {

    condvec <- unique(expdf$condition)
    resdflist <- parallel::mclapply(bytranslist, function(transtab, condvec,
        controlcondname, stresscondname) {

        ## Retrieve the column names for each comparison
        idxctrl <- grep(controlcondname, condvec)
        .checkempty(idxctrl, controlcondname)
        idxstress <- grep(stresscondname, condvec)
        .checkempty(idxstress, stresscondname)
        name1 <- paste0("mean_Fx_", condvec[idxctrl])  # nolint
        name2 <- paste0("mean_Fx_", condvec[idxstress])  # nolint
        diffname <- paste0("Diff_meanFx_", condvec[idxstress], "_",  # nolint
          condvec[idxctrl])

        ## Perform a kolmogorov-smirnoff test between the two columns
        resks <- suppressWarnings(stats::ks.test(transtab[, name1],
          transtab[, name2]))

        ## Calculate the area under the curve of the difference of means
        ## -> delta AUC
        deltadauc <- pracma::trapz(transtab[, "coord"], transtab[, diffname])
        ## Retrieve the p-value
        pvaldeltadaucks <- resks$p.value
        ## The KS test statistic is defined as the maximum value of the
        ## difference between A and B’s cumulative distribution functions (CDF)
        statdeltadaucks <- resks$statistic

        ## Build a one line data.frame with the proper col names
        ksdeltadaucdf <- data.frame(deltadauc, pvaldeltadaucks, statdeltadaucks)
        prefixvec <- c("dAUC", "p_dAUC", "D_dAUC")
        colnames(ksdeltadaucdf) <- paste(prefixvec, diffname, sep = "_")

        ## Retrieving transcript information
        infodf <- .returninfodf(transtab, nbwindows)

        ## Combining the two df as result
        resdf <- cbind(infodf, ksdeltadaucdf)
        return(resdf)
    }, condvec, controlcondname, stresscondname, mc.cores = nbcpu)

    resdf <- do.call("rbind", resdflist)

    ## Correct p-values using FDR
    idx <- grep("p_dAUC", colnames(resdf))
    fdrvec <- stats::p.adjust(resdf[, idx], method = "fdr")

    resdf <- cbind(resdf, fdrvec)
    colnamevec <- colnames(resdf)
    idxfdr <- which(colnamevec == "fdrvec")
    colnames(resdf)[idxfdr] <- paste0("adjFDR_", colnamevec[idx])  # nolint
    return(resdf)
}

.buildaucdf <- function(transtab, difffxname, resks, meanvalname,
  currentcond, nbwindows) {
    auc <- pracma::trapz(transtab[, "coord"], transtab[, difffxname])
    pvalaucks <- resks$p.value
    stataucks <- resks$statistic
    meanvaluefull <- mean(transtab[, meanvalname])
    aucdf <- data.frame(auc, pvalaucks, stataucks, meanvaluefull)
    prefixvec <- c("AUC", "p_AUC", "D_AUC", "MeanValueFull")
    colnames(aucdf) <- paste(prefixvec, currentcond, sep = "_")
    rownames(aucdf) <- NULL
    transinfo <- .returninfodf(transtab, nbwindows)
    res <- cbind(transinfo, aucdf)
    return(res)
}

.auc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1) {

  cumulative <- seq(1, nbwindows) / nbwindows
  condvec <- unique(expdf$condition)

  resdflist <- parallel::mclapply(bytranslist, function(transtab, condvec,
    cumulative, nbwindows) {
      ## Computing AUC, pval, and stat for each condition
      resauclist <- lapply(condvec, function(currentcond, transtab,
        cumulative, nbwindows) {
          ## Definition of column names
          difffxname <- paste0("diff_Fx_", currentcond) # nolint
          meanvalname <- paste0("mean_value_", currentcond) # nolint
          meanfxname <- paste0("mean_Fx_", currentcond) # nolint

          ## Perform a kolmogorov-smirnoff test between mean_Fx and cum.density
          resks <- suppressWarnings(stats::ks.test(transtab[, meanfxname],
            cumulative))
          ## Build data.frame with auc information for the current transcript
          aucdf <- .buildaucdf(transtab, difffxname, resks, meanvalname,
            currentcond, nbwindows)
          return(aucdf)
      }, transtab, cumulative, nbwindows)
      aucdf <- do.call("cbind", resauclist)
      return(aucdf)
  }, condvec, cumulative, nbwindows, mc.cores = nbcpu)

  aucallconditions <- do.call("rbind", resdflist)
  idxdup <- which(duplicated(colnames(aucallconditions)))
  if (!isTRUE(all.equal(length(idxdup), 0)))
    aucallconditions <- aucallconditions[, -idxdup]

  ## Correcting p-val with FDR
  idxpvalvec <- grep("p_AUC", colnames(aucallconditions))
  fdrlist <- lapply(idxpvalvec, function(idxpval, tab) {
    return(stats::p.adjust(tab[, idxpval], method = "fdr"))
  }, aucallconditions)
  fdrdf <- do.call("cbind", fdrlist)
  colnames(fdrdf) <- paste0("adjFDR_", colnames(aucallconditions)[idxpvalvec]) # nolint

  aucallconditions <- cbind(aucallconditions, fdrdf)
  return(aucallconditions)
}


#' Calculate Area Under Curve (AUC) and Differences of AUC for Transcript Data
#'
#' @description
#' This function computes the Area Under Curve (AUC) and the differences of AUC
#' between two conditions for a list of transcript data. It supports parallel
#' computation for efficiency. If only one condition is given, the differences
#' are not computed.
#'
#' @usage
#' allauc(bytranslistmean, expdf, nbwindows, nbcpu = 1,
#'  controlcondname = "ctrl", stresscondname = "HS", showtime = FALSE,
#'  verbose = TRUE)
#'
#' @param bytranslistmean A list of data frames, each containing transcript
#'                        level data with mean values for one or more
#'                        conditions.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param nbwindows An integer specifying the number of windows to consider for
#'                  AUC calculations.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'               parallel processing on bytranslistmean. Defaults to \code{1}.
#' @param controlcondname A string specifying the name of the control condition
#'                         Defaults to \code{"ctrl"}.
#' @param stresscondname A string specifying the name of the stress condition.
#'                        Defaults to \code{"HS"}.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical value indicating whether to print progress messages
#'                 Defaults to \code{TRUE}.
#'
#' @return A data frame containing the AUC and dAUC results for each transcript,
#'         along with associated statistical information.
#'
#' @details The function first checks if exactly two conditions are present in
#'          `expdf`. If so, it computes the differences in AUC between the two
#'          conditions using a Kolmogorov-Smirnov test. It then calculates the
#'          AUC for all conditions against a reference line (y=x). Results are
#'          merged by transcript and include adjusted p-values.
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#' 
#' ## Reading tables
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' 
#' ## Computing intermediate steps
#' avfilt <- averageandfilterexprs(expdf, transdf, expthres,
#'         showtime = FALSE, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' meandiff <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#' bytranslistmean <- split(meandiff, factor(meandiff$transcript))
#' 
#' ## Testing allauc
#' res <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
#'
#' @seealso
#' [genesECDF]
#'
#' @importFrom dplyr filter mutate select distinct
#' @importFrom parallel mclapply
#' @importFrom stats ks.test
#' @importFrom stats p.adjust
#' @importFrom magrittr %>%
#'
#' @export

allauc <- function(bytranslistmean, expdf, nbwindows, nbcpu = 1,
  controlcondname = "ctrl", stresscondname = "HS", showtime = FALSE,
  verbose = TRUE) {

    if (showtime) start_time <- Sys.time()
    if (verbose) message("\n\t ## Computing AUC and difference")

    if (isTRUE(all.equal(length(unique(expdf$condition)), 2))) {
        if (verbose) message("\t Computing the differences (d or delta) of AUC")
        daucallcond <- .dauc_allconditions(bytranslistmean, expdf, nbwindows,
          nbcpu, controlcondname, stresscondname)
    } else {
        if (verbose) message("\t dAUC not performed, only one condition",
          " submitted.")
    }

    ## Calculate the Area Under Curve (AUC), All conditions vs y=x
    ## Calculate Mean Value over the full gene body in All conditions.
    if (verbose) message("\t Computing the Area Under Curve (AUC)")
    aucallcond <- .auc_allconditions(bytranslistmean, expdf, nbwindows,
      nbcpu = nbcpu)

    ## Merging the two tables by transcript
    if (!isTRUE(all.equal(length(unique(expdf$condition)), 1))) {

      if (verbose) message("\t Merging results")
      allauc <- merge(aucallcond, daucallcond,
        by = c("gene", "transcript", "strand", "window_size"))
    } else {
      allauc <- aucallcond
    }

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }
    return(allauc)
}
