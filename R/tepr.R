#' Perform the tepr differential nascent rna-seq analysis
#'
#' @description
#' This function wraps the different steps of tepr to identify transcripts
#' with a significantly different nascent rna-seq signal.
#'
#' @usage
#' tepr(expdf, alldf, expthres, nbcpu = 1, rounding = 10, dontcompare = NULL,
#' controlcondname = "ctrl", stresscondname = "HS", replaceval = NA, pval = 0.1,
#' significant = FALSE, windsizethres = 50, countnathres = 20,
#' meanctrlthres = 0.5, meanstressthres = 0.5, pvaltheorythres = 0.1,
#' aucctrlthreshigher = -10, aucctrlthreslower = 15, aucstressthres = 15,
#' attenuatedpvalksthres = 2, outgrouppvalksthres = 0.2, showtime = FALSE,
#' verbose = TRUE)
#'
#' @param expdf A data frame containing experiment data that should have
#'  columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param alldf A data frame containing all transcript-related information,
#'  including biotype, chromosome, coordinates, transcript, gene, strand,
#'  window, ID and scores retrieved from the bedgraph files.
#' @param expthres A numeric value specifying the expression threshold.
#'  Transcripts with average expression values below this threshold will be
#'  filtered out from the returned transcript vector.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'  parallel computation on transcripts. The number of transcripts is equal to
#'  the number of lines provided as input of 'averageandfilterexprs'.
#'  Defaults to \code{1}.
#' @param rounding An integer specifying the rounding factor for computing ECDF.
#'  Default is \code{10}.
#' @param dontcompare An optional parameter to specify any conditions to exclude
#'  from the comparison. Defaults to \code{NULL}.
#' @param controlcondname A string specifying the name of the control condition
#'  Defaults to \code{"ctrl"}.
#' @param stresscondname A string specifying the name of the stress condition.
#'  Defaults to \code{"HS"}.
#' @param replaceval A value to replace non-significant attenuation values
#'  Defaults to \code{NA}.
#' @param pval A numeric value specifying the p-value threshold for significance
#'  of the KS test. Defaults to \code{0.1}.
#' @param significant A logical indicating whether to filter out non-significant
#'  attenuation values. Defaults to \code{FALSE}.
#' @param windsizethres A numeric threshold for the minimum window size. Default
#'  is 50.
#' @param countnathres A numeric threshold for the maximum number of missing
#'  data points for an experiment (NA values). Default is 20.
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
#'  processing should be indicated before ending. Defaults to \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'  Defaults to \code{TRUE}.
#'
#' @return
#' A list of data.frame made of two elements. The first data.frame is the
#' result of the function 'meandifference': For each condition, it provides
#' the mean values for the "value" and "Fx" columns (e.g. mean_value_ctrl, 
#' mean_Fx_ctrl columns) as the differences between the Fx column and coordinate
#' ratios (e.g., diff_Fx_ctrl column). The second data.frame is the result of
#' the function 'universegroup': It contains columns concerning AUC, knee
#' position, attenuation information, and columns defining the universe and
#' groups (see details).
#'






tepr <- function(expdf, alldf, expthres, nbcpu = 1, rounding = 10,
    dontcompare = NULL, controlcondname = "ctrl", stresscondname = "HS",
    replaceval = NA, pval = 0.1, significant = FALSE, windsizethres = 50,
    countnathres = 20, meanctrlthres = 0.5, meanstressthres = 0.5,
    pvaltheorythres = 0.1, aucctrlthreshigher = -10, aucctrlthreslower = 15,
    aucstressthres = 15, attenuatedpvalksthres = 2, outgrouppvalksthres = 0.2,
    showtime = FALSE, verbose = TRUE) {

    ## This function calculates the average expression levels for transcripts
    ## from a provided expression data frame and filters out transcripts based
    ## on a specified expression threshold.
    resallexprs <- averageandfilterexprs(expdf, alldf, expthres, showtime,
        verbose = TRUE)

    ## This function takes a list of expression data frames, a condition
    ## information data frame, and counts the number of NA values for each
    # transcript based on strand and condition.
    rescountna <- countna(resallexprs, expdf, nbcpu, showtime, verbose)

    ## This function calculates the empirical cumulative distribution function
    ## (ECDF) for expressed genes across multiple transcripts.
    resecdflist <- genesECDF(resallexprs, expdf, nbcpu, rounding, showtime,
        verbose)
    resecdf <- resecdflist[[1]]
    nbwindows <- resecdflist[[2]]

    ## This function calculates the mean values, mean Fx (ECDF) and ECDF
    ## differences (Fx) for expression data, across different experimental
    ## conditions.
    resmeandiff <- meandifference(resecdf, expdf, nbwindows, showtime, verbose)

    ## Split the results by transcripts
    if (showtime) start_time_split <- Sys.time()
    if (verbose) message("Split the results by transcripts")
    bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
    if (showtime) {
        end_time_split <- Sys.time()
        timing <- end_time_split - start_time_split
        message("\t\t ## Analysis performed in: ", format(timing, digits = 2))
    }

    ## This function computes the Area Under Curve (AUC) and the differences of
    ## AUC between two conditions for a list of transcript data.
    resauc <- allauc(bytranslistmean, expdf, nbwindows, nbcpu,
        dontcompare, controlcondname, stresscondname, showtime, verbose)

    ## This function identifies the knee point (i.e., point of maximum change)
    ## and the maximum difference in the empirical cumulative distribution
    ## function (ECDF) for each transcript, across different experimental
    ## conditions.
    resknee <- kneeid(bytranslistmean, expdf, nbcpu, showtime, verbose)

    ## This function computes the attenuation values for each window of each
    ## transcript based on the data frames obtained with the functions 'allauc',
    ## 'kneeid', and 'countna'.
    resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
        resmeandiff, nbcpu, significant, replaceval, pval, showtime, verbose)

    ## This function categorizes genes into a "Universe" and assigns them into
    ## groups such as "Attenuated" or "Outgroup" based on transcription data and
    ## thresholds.
    res <- universegroup(resatt, controlcondname, stresscondname, windsizethres,
        countnathres, meanctrlthres, meanstressthres, pvaltheorythres,
        aucctrlthreshigher, aucctrlthreslower, aucstressthres,
        attenuatedpvalksthres, outgrouppvalksthres, showtime, verbose)

    ## Return variables necessary for plotting
    return(list(resmeandiff, res))
}