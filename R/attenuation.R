.summarytrans <- function(bytransmeanlist, nbcpu) {
  summarydflist <- parallel::mclapply(bytransmeanlist, function(trans) {
    coor1 <- min(trans$coor1)
    coor2 <- max(trans$coor2)
    return(data.frame(chr = trans$chr[1], coor1, coor2,
          strand = trans$strand[1], gene = trans$gene[1],
          transcript = trans$transcript[1], size = coor2 - coor1 + 1))
  }, mc.cores = nbcpu)
  summarydf <- do.call("rbind", summarydflist)
  return(summarydf)
}

.computeupdown <- function(completbytrans, condvec, nbcpu) {

  updownbytranslist <- parallel::mclapply(completbytrans,
    function(trans, condvec) {

    ## Ordering by coordinates (security)
    trans <- trans[order(trans$coord), ]

    ## For each condition
    updownlist <- lapply(condvec, function(cond, trans) {
      kneecolname <- paste0("knee_AUC_", cond)
      meancolname <- paste0("mean_value_", cond)

      idxup <- which(trans$coord <= trans[, kneecolname])
      if (isTRUE(all.equal(length(idxup), 0)))
        stop("\n\t Problem in retrieving idxup, this should not happen. ",
            "Contact the developer.\n")
      upmean <- mean(trans[idxup, meancolname])

      idxdown <- which(trans$coord >= trans[, kneecolname] &
                          trans$coord <= max(trans$coord))
      if (isTRUE(all.equal(length(idxdown), 0)))
        stop("\n\t Problem in retrieving idxdown, this should not happen. ",
            "Contact the developer.\n")
      downmean <- mean(trans[idxdown, meancolname])

      ## Calculating attenuation
      att <- 100 - downmean / upmean * 100

      res <- data.frame(trans$transcript[1], upmean, downmean, att)
      colnames(res) <- c("transcript", paste0("UP_mean_", cond),
              paste0("DOWN_mean_", cond), paste0("Attenuation_", cond))
      return(res)
    }, trans)

    return(do.call("cbind", updownlist))
  }, condvec, mc.cores = nbcpu)

  updowndf <- do.call("rbind", updownbytranslist)
  idxdup <- which(duplicated(colnames(updowndf)))
  if (!isTRUE(all.equal(length(idxdup), 0)))
    updowndf <- updowndf[, -idxdup]

  return(updowndf)
}

.filterattenuation <- function(auckneenasumatt, condvec, pval, replaceval, # nolint
    verbose) {

        mat <- auckneenasumatt
        if (verbose) message("\t\t Replacing non-significant attenuations by ",
            replaceval)
        invisible(sapply(condvec, function(cond, replaceval) {
            pauccond <- paste0("p_AUC_", cond)
            ## Replacing Attenuation value if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("Attenuation_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("Attenuation_", cond)]]))
            ## Replacing knee values if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("knee_AUC_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("knee_AUC_", cond)]]))
        }, replaceval))

        return(mat)
}

#' Calculate Attenuation from AUC and Other Transcript Features
#'
#' @description
#' This function computes the attenuation values for each window of each
#' transcript based on the data frames obtained with the functions 'allauc',
#' 'kneeid', and 'countna'.
#'
#' @usage
#' attenuation(allaucdf, kneedf, matnatrans, bytranslistmean, expdf, dfmeandiff,
#' nbcpu = 1, significant = FALSE, replaceval = NA, pval = 0.1,
#' showtime = FALSE, verbose = TRUE)
#'
#' @param allaucdf A data frame containing AUC results for transcripts (see
#'                 allauc).
#' @param kneedf A data frame containing the inflection points (see kneeid).
#' @param matnatrans A data frame containing the number of missing values per
#'                   transcript (see countna).
#' @param bytranslistmean A list of data frames with mean values by transcripts.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param dfmeandiff A data frame containing means and differences in mean
#'                  values, if more than one condition. (see meandifference).
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'               parallel processing. The parallelization is done on
#'               bytranslistmean whose number of elements is equal to the
#'               number of lines provided as input of 'averageandfilterexprs'.
#'               Defaults to \code{1}.
#' @param significant A logical indicating whether to filter out non-significant
#'                    attenuation values. Defaults to \code{FALSE}.
#' @param replaceval A value to replace non-significant attenuation values
#'                   Defaults to \code{NA}.
#' @param pval A numeric value specifying the p-value threshold for significance
#'              of the KS test. Defaults to \code{0.1}.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical value indicating whether to print progress messages
#'                 Defaults to \code{TRUE}.
#'
#' @return A data frame containing the computed attenuation values along with
#'         associated transcript information.
#'
#' @details The function merges several data frames to create a comprehensive
#'          dataset for each transcript. It computes mean values for the "up"
#'          and "down" segments of the transcript. The direction is determined
#'          by comparing the coordinates to the knee values. up = coord < knee
#'          and down = coord > knee. The up and down indexes are then retrieved
#'          and the attenuation scores are computed as:
#'              att <- 100 - downmean / upmean * 100
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
#'        showtime = FALSE, verbose = FALSE)
#' rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' resmeandiff <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#' bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
#' resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
#' resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
#' 
#' ## Testing attenuation
#' resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
#'         resmeandiff, verbose = FALSE)
#'
#' @seealso
#' [allauc()], [kneeid()], [countna()], [meandifference()]
#'
#' @importFrom dplyr filter mutate select
#' @importFrom parallel mclapply
#' @importFrom stats ks.test
#' @importFrom stats p.adjust
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#'
#' @export

attenuation <- function(allaucdf, kneedf, matnatrans, bytranslistmean, expdf,
  dfmeandiff, nbcpu = 1, significant = FALSE, replaceval = NA, pval = 0.1,
  showtime = FALSE, verbose = TRUE) {

      if (showtime) start_time <- Sys.time()
      if (verbose) message("\n\t ## Calculating attenuation")
      if (verbose) message("\t Merging tables")
      allaucknee <- merge(allaucdf, kneedf, by = "transcript")
      mergecolnames <- c("gene", "transcript", "strand")
      allauckneena <- merge(allaucknee, matnatrans, by = mergecolnames)

      if (verbose) message("\t Building summary")
      summarydf <- .summarytrans(bytranslistmean, nbcpu)
      if (verbose) message("\t Merging summary")
      auckneenasum <- merge(summarydf, allauckneena, by = mergecolnames)

      ## Merging the mean table with the previous one
      if (verbose) message("\t Merging detailed mean table with summary")
      complet <- merge(dfmeandiff, auckneenasum, by = mergecolnames)

      ## Splitting the previous table by transcript
      if (verbose) message("\t Splitting the previous table by transcript")
      completbytrans <- split(complet, factor(complet$transcript))
      condvec <- unique(expdf$condition)

      ## For each transcript
      if (verbose) message("\t Computing up and down mean")
      updowndf <- .computeupdown(completbytrans, condvec, nbcpu)

      ## Merging attenuation to the complete table
      if (verbose) message("\t Merging attenuation to the complete table")
      auckneenasumatt <- merge(auckneenasum, updowndf, by = "transcript")

      ## Replace the attenuation values by replaceval if p_AUC_cond >= pval
      if (significant) {
        if (verbose) message("\t Keeping significant attenuation")
        auckneenasumatt <- .filterattenuation(auckneenasumatt, condvec, pval,
            replaceval, verbose)
      }

      if (showtime) {
        end_time <- Sys.time()
        timing <- end_time - start_time
        message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
      }

      return(auckneenasumatt)
}
