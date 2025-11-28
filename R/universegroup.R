#' Define Universe and Group of Genes Based on Expression Data
#'
#' @description
#' This function categorizes genes into a "Universe" and assigns them into
#' groups such as "Attenuated" or "Outgroup" based on transcription data and
#' thresholds. The universe is defined by thresholds for window size, missing
#' data count, mean transcription levels, and p-values. Genes are further
#' classified into groups based on conditions related to AUC and p-value
#' thresholds.
#'
#' @usage
#' universegroup(completedf, expdf, controlname = "ctrl", stressname = "HS",
#' windsizethres = 50, countnathres = 20, meanctrlthres = 0.5,
#' meanstressthres = 0.5, pvaltheorythres = 0.1, aucctrlthreshigher = -10,
#' aucctrlthreslower = 15, aucstressthres = 15, attenuatedpvalksthres = 2,
#' outgrouppvalksthres = 0.2, showtime = FALSE, verbose = TRUE)
#'
#' @param completedf A data frame obtained with the function attenuation.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param controlname A string representing the control condition name. Default
#'  is \code{"ctrl"}.
#' @param stressname A string representing the stress condition name. Default
#'  is \code{"HS"}.
#' @param windsizethres A numeric threshold for the minimum window size. Default
#'  is 50.
#' @param countnathres A numeric threshold for the maximum number of missing
#'  data points (NA values). Default is 20.
#' @param meanctrlthres A numeric threshold for the minimum mean transcription
#'  value in the control condition. Default is 0.5.
#' @param meanstressthres A numeric threshold for the minimum mean transcription
#'  value in the stress condition. Default is 0.5.
#' @param pvaltheorythres A numeric threshold for the minimum p-value used to
#'  define the universe of genes. Default is 0.1.
#' @param aucctrlthreshigher A numeric threshold for the lower bound of the
#'  control AUC value in the outgroup classification. Default is -10.
#' @param aucctrlthreslower A numeric threshold for the upper bound of the
#'  control AUC value in the outgroup classification. Default is 15.
#' @param aucstressthres A numeric threshold for the minimum stress AUC value
#'  used to classify attenuated genes. Default is 15.
#' @param attenuatedpvalksthres A numeric threshold for the negative log10 of
#'  the p-value (from KS test) for defining attenuated genes. Default is 2.
#' @param outgrouppvalksthres A numeric threshold for the maximum KS p-value
#'  used to define the outgroup. Default is 0.2.
#' @param showtime A logical value indicating if the duration of the function
#'                  processing should be indicated before ending. Defaults to
#'                  \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return A modified data frame with two additional columns: \code{Universe},
#'  indicating whether each gene is part of the universe, and \code{Group},
#'  classifying the genes into groups such as "Attenuated", "Outgroup", or
#'  \code{NA}.
#'
#' @details
#' A transcript belongs to "Universe" if:
#'      window_size > windsizethres & Count_NA < countnathres &
#'      meanctrl > meanctrlthres & meanstress > meanstressthres &
#'       pvaltheory > pvaltheorythres
#'
#' If only one condition is provided, a transcript belongs to "Universe" if:
#'      window_size > windsizethres & Count_NA < countnathres &
#'      meanctrl > meanctrlthres & pvaltheory > pvaltheorythres
#'
#' A transcript belongs to the groups:
#' - \strong{Attenuated}: if Universe == TRUE & aucstress > aucstressthres & -log10(pvalks) > attenuatedpvalksthres
#' - \strong{Outgroup}: if Universe == TRUE & pvalks > outgrouppvalksthres & aucctrl > aucctrlthreshigher & aucctrl < aucctrlthreslower
#'
#' If only one condition is provided:
#' - \strong{Attenuated}: if Universe == TRUE & aucctrl > aucctrlthreslower
#' - \strong{Outgroup}: if Universe == TRUE & aucctrl > aucctrlthreshigher & aucctrl < aucctrlthreslower
#'
#' This function is useful for classifying genes in transcriptomics data based
#' on their transcriptional response to different experimental conditions.
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
#'        showtime = FALSE, verbose = FALSE)
#' rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' resmeandiff <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#' bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
#' resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
#' resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
#' resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
#'         resmeandiff, verbose = FALSE)
#' ## Testing universegroup
#' resug <- universegroup(resatt, expdf, verbose = FALSE)
#'
#' @seealso
#' [attenuation]
#'
#' @importFrom dplyr mutate relocate select filter
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#'
#' @export

universegroup <- function(completedf, expdf, controlname = "ctrl", # nolint
    stressname = "HS", windsizethres = 50, countnathres = 20,
    meanctrlthres = 0.5, meanstressthres = 0.5, pvaltheorythres = 0.1,
    aucctrlthreshigher = -10, aucctrlthreslower = 15, aucstressthres = 15,
    attenuatedpvalksthres = 2, outgrouppvalksthres = 0.2, showtime = FALSE,
    verbose = TRUE) {

    if (showtime) start_time <- Sys.time()
    if (verbose) message("\n\t ## Computing universe and group columns")

    meanctrl <- paste("MeanValueFull", controlname, sep = "_")
    meanstress <- paste("MeanValueFull", stressname, sep = "_")
    pvaltheory <- paste("adjFDR_p_AUC", controlname, sep = "_")
    aucctrl <- paste("AUC", controlname, sep = "_")
    aucstress <- paste("AUC", stressname, sep = "_")
    pvalks <- paste0("adjFDR_p_dAUC_Diff_meanFx_", stressname, "_", controlname)
    condvec <- unique(expdf$condition)

    ## Computing the Universe column: If only one condition is provided, only
    ## control columns are used
    if (!isTRUE(all.equal(length(condvec), 1))) {
        completedf <- completedf %>%
        dplyr::mutate(Universe = ifelse(
            .data$window_size > windsizethres &
            .data$Count_NA < countnathres &
            !!sym(meanctrl) > meanctrlthres & # nolint
            !!sym(meanstress) > meanstressthres &
            !!sym(pvaltheory) > pvaltheorythres, TRUE, FALSE)) %>%
            dplyr::relocate("Universe", .before = 1)  # nolint
    } else {
        completedf <- completedf %>%
        dplyr::mutate(Universe = ifelse(
            .data$window_size > windsizethres &
            .data$Count_NA < countnathres &
            !!sym(meanctrl) > meanctrlthres & # nolint
            !!sym(pvaltheory) > pvaltheorythres, TRUE, FALSE)) %>%
            dplyr::relocate("Universe", .before = 1)  # nolint
    }

    ## Computing the Group column
    if (!isTRUE(all.equal(length(condvec), 1))) {
        completedf <- completedf %>%
        dplyr::mutate(
            Group = ifelse(.data$Universe == TRUE &
                !!sym(aucstress) > aucstressthres &
                -log10(!!sym(pvalks)) > attenuatedpvalksthres, "Attenuated",
                NA),
            Group = ifelse(.data$Universe == TRUE &
                !!sym(pvalks) > outgrouppvalksthres &
                !!sym(aucctrl) > aucctrlthreshigher &
                !!sym(aucctrl) < aucctrlthreslower, "Outgroup",
                    .data$Group)) %>%
                dplyr::relocate("Group", .before = 2)
    } else {
        completedf <- completedf %>%
        dplyr::mutate(
            Group = ifelse(.data$Universe == TRUE &
                !!sym(aucctrl) > aucctrlthreslower, "Attenuated",
                NA),
            Group = ifelse(.data$Universe == TRUE &
                !!sym(aucctrl) > aucctrlthreshigher &
                !!sym(aucctrl) < aucctrlthreslower, "Outgroup",
                    .data$Group)) %>% dplyr::relocate("Group", .before = 2)
    }

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }

    return(completedf)
}
