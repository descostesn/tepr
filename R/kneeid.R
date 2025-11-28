.retrievekneeandmax <- function(condvec, transtable) { # nolint

reslist <- lapply(condvec, function(cond, transtable) {
        difffxname <- paste0("diff_Fx_", cond)
        difffxvec <- transtable[, difffxname]
        ## If equality of difference within the same gene it takes the closest
        ## knee from the TSS # nolint
        resrow <- transtable[which(difffxvec == max(difffxvec)), ] %>% # nolint
            dplyr::slice_min(.data$coord, n = 1) # nolint
        res <- data.frame(resrow$coord, resrow[, difffxname])
        colnames(res) <- c(paste0("knee_AUC_", cond), paste0("max_",
            difffxname))
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
#' @param expdf A data frame containing experiment data that should have
#'   columns named 'condition', 'replicate', 'strand', and 'path'.
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
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Calculating necessary results
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' avfilt <- averageandfilterexprs(expdf, transdf, expthres,
#'         showtime = FALSE, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' meandiff <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#' bytranslistmean <- split(meandiff, factor(meandiff$transcript))
#'
#' ## Testing kneeid
#' reskneeid <- kneeid(bytranslistmean, expdf, verbose = FALSE)
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
  if (verbose) message("\n\t ## Calculating knee")
  condvec <- unique(expdf$condition)
  bytransres <- parallel::mclapply(transdflist, function(transtable, condvec) {
      bycondreslist <- .retrievekneeandmax(condvec, transtable)
       return(cbind(transcript = transtable$transcript[1],
        do.call("cbind", bycondreslist)))
    }, condvec, mc.cores = nbcpu)
  res <- do.call("rbind", bytransres)

  if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
  }

  return(res)
}




.alldfcond <- function(currentexpdf, alldf) {

    ## Building vectors with the column names specific to the two conditions
    namecols <- paste0(currentexpdf$condition, "_rep", currentexpdf$replicate,
            ".", currentexpdf$strand)
    idxcolcond <- unlist(lapply(namecols,
            function(x, alldf) grep(x, colnames(alldf)), alldf))

    ## Limiting alldf to the two defined conditions
    alldfcond <- alldf[, c(seq_len(9), idxcolcond)]
    return(alldfcond)
}

#' Calculate knee for each condition separately
#'
#' @description
#' This function calculate knees for each condition separately
#'
#' @usage
#' kneeallconds(alldf, expdf, expthres, nbcpu = 1, rounding = 10,
#' showtime = FALSE, verbose = TRUE)
#'
#' @param alldf A data frame containing all transcript-related information,
#'  including biotype, chromosome, coordinates, transcript, gene, strand,
#'  window, ID and scores retrieved from the bedgraph files.
#' @param expdf A data frame containing experiment data that should have
#'  columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param expthres A numeric value specifying the expression threshold.
#'  Transcripts with average expression values below this threshold will be
#'  filtered out from the returned transcript vector.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'  parallel computation on transcripts. The number of transcripts is equal to
#'  the number of lines provided as input of 'averageandfilterexprs'.
#'  Defaults to \code{1}.
#' @param rounding An integer specifying the rounding factor for computing ECDF.
#'  Default is \code{10}.
#' @param showtime A logical value indicating if the duration of the function
#'  processing should be indicated before ending. Defaults to \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return
#' A data.frame with the columns \code{transcript}, \code{gene}, and
#' \code{strand window_size}. For each condition 'cond': \code{AUC_cond},
#' \code{p_AUC_cond}, \code{D_AUC_cond}, \code{MeanValueFull_cond},
#' \code{adjFDR_p_AUC_cond}, \code{knee_AUC_cond}, and \code{max_diff_Fx_cond}.
#'
#' @details
#' The kneeallconds function calls successively for each condition:
#' \itemize{
#'   \item averageandfilterexprs This function calculates the average expression levels for transcripts from alldf that was obtained with the 'preprocessing' function. It filters out transcripts based on the 'expthres' expression threshold. The function also renames the columns in the output data frame to include mean expression values. It returns a list containing the original alldf with the mean columns added and a character vector of transcripts that meet the filtering criteria.
#'   \item genesECDF It takes the result of 'averageandfilterexprs' as input and a) filters the main expression table to retain only the expressed transcripts; b) Splits the data by each transcript; c) For each transcript, computes ECDF values for the score columns while respecting the strand orientation ("plus" or "minus"); d) Combines the ECDF results into a final data frame. It returns a list containing two elements: A data frame with ECDF results for each transcript (concatdf) and an integer indicating the number of rows in each transcript table.
#'   \item meandifference It takes the result of 'genesECDF' as input and calculates the mean values, mean Fx (ECDF) and ECDF differences (Fx) for expression data, across different experimental conditions. It returns a data frame that contains, for each condition: mean values for the "value" and "Fx" columns (e.g., \code{mean_value_ctrl}, \code{mean_Fx_ctrl}) and the differences between the \code{Fx} column and coordinate ratios (e.g., \code{diff_Fx_ctrl}).
#'   \item Splitting Split the results of 'meandifference' by transcripts and stores the list into bytranslistmean.
#'   \item allauc It uses the previously computed variable 'bytranslistmean' and it computes the Area Under Curve (AUC) and the differences of AUC between two conditions for a list of transcript data. It returns a data frame containing the AUC and dAUC results for each transcript, along with associated statistical information.
#'   \item kneeid It uses the previously computed variable 'bytranslistmean' and identifies the knee point (i.e., point of maximum change) and the maximum difference in the empirical cumulative distribution function (ECDF) for each transcript, across different experimental conditions. It returns a data frame where each row corresponds to a transcript and contains the coordinates of the knee point and the maximum ECDF difference for each condition.
#' }
#'
#' @examples
#' exptabpath <- system.file("extdata", "exptab.csv", package="tepr")
#' alldfpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expdf <- read.csv(exptabpath)
#' alldf <- read.delim(alldfpath, header = FALSE)
#' expthres <- 0.1
#' kneedf <- kneeallconds(alldf, expdf, expthres, verbose = FALSE)
#'
#' @seealso
#' [averageandfilterexprs()], [genesECDF()], [meandifference()],
#' [kneeid()]
#'
#' @importFrom purrr reduce
#'
#' @export

kneeallconds <- function(alldf, expdf, expthres, nbcpu = 1, rounding = 10,
    showtime = FALSE, verbose = TRUE) {

    if (showtime) start_knee <- Sys.time()

    checkexptab(expdf)

    ## Building col names for alldf
    alldf <- .buildcolnames(expdf, alldf)

    ## Splitting expdf by conditions
    expdfcondlist <- split(expdf, factor(expdf$condition))

    ## Calling building of knee for each comparison of matcond
    kneelist <- lapply(expdfcondlist, function(currentexpdf, alldf, expthres,
        nbcpu, rounding, showtime, verbose) {

            if (verbose) message("\n ## Calculating knee for ",
                unique(currentexpdf$condition))
            alldfcond <- .alldfcond(currentexpdf, alldf)

            resallexprs <- averageandfilterexprs(currentexpdf, alldfcond,
                expthres, showtime, verbose)
            resecdflist <- genesECDF(resallexprs, nbcpu, rounding,
                showtime, verbose)
            resmeandiff <- meandifference(resecdflist[[1]], currentexpdf,
                resecdflist[[2]], showtime, verbose)
            if (verbose) message("\n\t ## Splitting by transcript")
            bytranslist <- split(resmeandiff, factor(resmeandiff$transcript))
            resauc <- allauc(bytranslist, currentexpdf,
                nbwindows = resecdflist[[2]], nbcpu = nbcpu,
                showtime = showtime, verbose = verbose)
            resknee <- kneeid(bytranslist, currentexpdf, nbcpu, showtime,
                verbose)
            resmerge <- merge(resauc, resknee, by = "transcript")

            rm(alldfcond, resallexprs, resecdflist, resmeandiff, bytranslist,
                resauc, resknee)
            invisible(gc())
            return(resmerge)

    }, alldf, expthres, nbcpu, rounding, showtime, verbose)

    ## Combine results
    if (verbose) message("\n Merging results into a single table.\n")
    joincols <- c("transcript", "gene", "strand", "window_size")
    kneedf <- purrr::reduce(kneelist, dplyr::left_join, by = joincols)

    if (showtime) {
      end_knee <- Sys.time()
      timing <- end_knee - start_knee
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }

    return(kneedf)
}
