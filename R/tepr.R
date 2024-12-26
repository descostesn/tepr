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