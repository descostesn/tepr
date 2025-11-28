#' Perform the tepr differential nascent rna-seq analysis
#'
#' @description
#' This function wraps the different steps of tepr to identify transcripts
#' with a significantly different nascent rna-seq signal.
#'
#' @usage
#' tepr(expdf, alldf, expthres, nbcpu = 1, rounding = 10,
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
#' @details
#' The tepr function calls successively:
#' \itemize{
#'   \item averageandfilterexprs This function calculates the average expression levels for transcripts from alldf that was obtained with the 'preprocessing' function. It filters out transcripts based on the 'expthres' expression threshold. The function also renames the columns in the output data frame to include mean expression values. It returns a list containing the original alldf with the mean columns added and a character vector of transcripts that meet the filtering criteria.
#'   \item countna It takes the result of 'averageandfilterexprs' as input and counts the number of NA values for each transcript based on strand and condition. NA represent missing scores that were filtered out from the black list and mappability track. It returns a data frame where each row corresponds to a transcript, along with its associated gene, strand, and the count of NA values.
#'   \item genesECDF It takes the result of 'averageandfilterexprs' as input and a) filters the main expression table to retain only the expressed transcripts; b) Splits the data by each transcript; c) For each transcript, computes ECDF values for the score columns while respecting the strand orientation ("plus" or "minus"); d) Combines the ECDF results into a final data frame. It returns a list containing two elements: A data frame with ECDF results for each transcript (concatdf) and an integer indicating the number of rows in each transcript table.
#'   \item meandifference It takes the result of 'genesECDF' as input and calculates the mean values, mean Fx (ECDF) and ECDF differences (Fx) for expression data, across different experimental conditions. It returns a data frame that contains, for each condition: mean values for the "value" and "Fx" columns (e.g., \code{mean_value_ctrl}, \code{mean_Fx_ctrl}) and the differences between the \code{Fx} column and coordinate ratios (e.g., \code{diff_Fx_ctrl}).
#'   \item Splitting Split the results of 'meandifference' by transcripts and stores the list into bytranslistmean.
#'   \item allauc It uses the previously computed variable 'bytranslistmean' and it computes the Area Under Curve (AUC) and the differences of AUC between two conditions for a list of transcript data. It returns a data frame containing the AUC and dAUC results for each transcript, along with associated statistical information.
#'   \item kneeid It uses the previously computed variable 'bytranslistmean' and identifies the knee point (i.e., point of maximum change) and the maximum difference in the empirical cumulative distribution function (ECDF) for each transcript, across different experimental conditions. It returns a data frame where each row corresponds to a transcript and contains the coordinates of the knee point and the maximum ECDF difference for each condition.
#'   \item attenuation It uses the results of the previous functions and it computes the attenuation values for each window of each transcript. It returns a data frame containing the computed attenuation values along with associated transcript information.
#'   \item universegroup Using the table produced by 'attenuation', it categorizes genes into a "Universe" and assigns them into groups such as "Attenuated" or "Outgroup" based on transcription data and thresholds. The universe is defined by thresholds for window size, missing data count, mean transcription levels, and p-values. Genes are further classified into groups based on conditions related to AUC and p-value thresholds. A transcript belongs to "Universe" if (window_size > windsizethres & Count_NA < countnathres & meanctrl > meanctrlthres & meanstress > meanstressthres & pvaltheory > pvaltheorythres). A transcript belongs to the groups: - \strong{Attenuated}: if Universe == TRUE & aucstress > aucstressthres & -log10(pvalks) > attenuatedpvalksthres. - \strong{Outgroup}: if Universe == TRUE & pvalks > outgrouppvalksthres & aucctrl > aucctrlthreshigher & aucctrl < aucctrlthreslower.
#' }
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Reading files
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#'
#' ## Testing tepr
#' res <- tepr(expdf, transdf, expthres, verbose = FALSE)
#'
#' @seealso
#' [averageandfilterexprs()], [countna()], [genesECDF()], [meandifference()],
#' [allauc()], [kneeid()], [attenuation()], [universegroup()], [preprocessing()]
#'
#' @export

tepr <- function(expdf, alldf, expthres, nbcpu = 1, rounding = 10,
    controlcondname = "ctrl", stresscondname = "HS",
    replaceval = NA, pval = 0.1, significant = FALSE, windsizethres = 50,
    countnathres = 20, meanctrlthres = 0.5, meanstressthres = 0.5,
    pvaltheorythres = 0.1, aucctrlthreshigher = -10, aucctrlthreslower = 15,
    aucstressthres = 15, attenuatedpvalksthres = 2, outgrouppvalksthres = 0.2,
    showtime = FALSE, verbose = TRUE) {

    if (showtime) start_tepr <- Sys.time()

    if (length(unique(expdf$condition)) > 2)
        stop("\n[tepr] Error: Too many conditions.\n",
            "  Found more than 2 conditions. Use teprmulti() instead.\n")

    ## This function calculates the average expression levels for transcripts
    ## from a provided expression data frame and filters out transcripts based
    ## on a specified expression threshold.
    resallexprs <- averageandfilterexprs(expdf, alldf, expthres, showtime,
        verbose)

    ## This function takes a list of expression data frames, a condition
    ## information data frame, and counts the number of NA values for each
    # transcript based on strand and condition.
    rescountna <- countna(resallexprs, expdf, nbcpu, showtime, verbose)

    ## This function calculates the empirical cumulative distribution function
    ## (ECDF) for expressed genes across multiple transcripts.
    resecdflist <- genesECDF(resallexprs, nbcpu, rounding, showtime,
        verbose)
    resecdf <- resecdflist[[1]]
    nbwindows <- resecdflist[[2]]

    ## This function calculates the mean values, mean Fx (ECDF) and ECDF
    ## differences (Fx) for expression data, across different experimental
    ## conditions.
    resmeandiff <- meandifference(resecdf, expdf, nbwindows, showtime, verbose)

    ## Split the results by transcripts
    if (showtime) start_time_split <- Sys.time()
    if (verbose) message("\n\t ## Split the results by transcripts")
    bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
    if (showtime) {
        end_time_split <- Sys.time()
        timing <- end_time_split - start_time_split
        message("\t\t -- Analysis performed in: ", format(timing, digits = 2))
    }

    ## This function computes the Area Under Curve (AUC) and the differences of
    ## AUC between two conditions for a list of transcript data.
    resauc <- allauc(bytranslistmean, expdf, nbwindows, nbcpu, controlcondname,
        stresscondname, showtime, verbose)

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
    res <- universegroup(resatt, expdf, controlcondname, stresscondname,
        windsizethres, countnathres, meanctrlthres, meanstressthres,
        pvaltheorythres, aucctrlthreshigher, aucctrlthreslower, aucstressthres,
        attenuatedpvalksthres, outgrouppvalksthres, showtime, verbose)

    if (showtime) {
        end_tepr <- Sys.time()
        timing <- end_tepr - start_tepr
        message("\n\t\t -- tepr analysis performed in: ",
            format(timing, digits = 2))
    }
    ## Return variables necessary for plotting
    reslist <- list(resmeandiff, res)
    return(reslist)
}





.restepr <- function(saveobjectpath, compname, reload, expdftwocond, alldftwocond,
    expthres, nbcpu, rounding, condonename, condtwoname, replaceval,
    pval, significant, windsizethres, countnathres, meancondonethres,
    meancondtwothres, pvaltheorythres, auccondonethreshigher, auccondonethreslower,
    auccondtwothres, attenuatedpvalksthres, outgrouppvalksthres, showtime,
    verbose) {

    filepathname <- file.path(saveobjectpath, paste0(compname, ".rds"))

    if (!reload || !file.exists(filepathname)) {

        restepr <- tepr(expdf = expdftwocond, alldf = alldftwocond,
            expthres = expthres, nbcpu = nbcpu, rounding = rounding,
            controlcondname = condonename, stresscondname = condtwoname,
            replaceval = replaceval, pval = pval, significant = significant,
            windsizethres = windsizethres, countnathres = countnathres,
            meanctrlthres = meancondonethres, meanstressthres = meancondtwothres,
            pvaltheorythres = pvaltheorythres,
            aucctrlthreshigher = auccondonethreshigher,
            aucctrlthreslower = auccondonethreslower,
            aucstressthres = auccondtwothres,
            attenuatedpvalksthres = attenuatedpvalksthres,
            outgrouppvalksthres = outgrouppvalksthres, showtime = showtime,
            verbose = verbose)

        names(restepr) <- c(paste("resmeandiff", compname, sep = "_"),
            paste("resunigroupatt", compname, sep = "_"))

        if (!is.na(saveobjectpath)) {
            if (verbose) message("\t\t Saving to ", filepathname)
            saveRDS(restepr, file = filepathname)
        }

    } else {
        if (verbose) message("\t\t\t Loading ", filepathname)
        restepr <- readRDS(filepathname)
    }

    return(restepr)
}

.reslist <- function(filepathname, matcond, verbose, expdf, alldf, expthres,
    nbcpu, rounding, replaceval, pval, significant, windsizethres,
    countnathres, meancondonethres, meancondtwothres, pvaltheorythres,
    auccondonethreshigher, auccondonethreslower, auccondtwothres,
    attenuatedpvalksthres, outgrouppvalksthres, saveobjectpath, reload,
    showtime, showmemory) {

        reslist <- apply(matcond, 2, function(currentcol, verbose, expdf, alldf,
        expthres, nbcpu, rounding, replaceval, pval, significant,
        windsizethres, countnathres, meancondonethres, meancondtwothres,
        pvaltheorythres, auccondonethreshigher, auccondonethreslower, auccondtwothres,
        attenuatedpvalksthres, outgrouppvalksthres, saveobjectpath,
        reload, showtime, showmemory) {

        condonename <- currentcol[1]
        condtwoname <- currentcol[2]
        compname <- paste(condonename, condtwoname, sep = "_vs_")
        if (verbose) message("\n\n Comparison of ", compname)

        ## Limiting expdf on the two defined conditions
        idxexp <- as.vector(sapply(currentcol, function(condname, expdf) {
            return(which(expdf$condition == condname))}, expdf))
        expdftwocond <- expdf[idxexp, ]
        rownames(expdftwocond) <- NULL

        ## Building vectors with the column names specific to the two conditions
        namecols <- paste0(expdftwocond$condition, "_rep", expdftwocond$replicate,
            ".", expdftwocond$strand)
        idxcoltwoconds <- unlist(lapply(namecols,
            function(x, alldf) grep(x, colnames(alldf)), alldf))

        ## The info columns are biotype, chr, coor1, coor2, transcript, gene,
        ## strand, window, id. This is reflected by seq_len(9)
        ## Limiting alldf to the two defined conditions
        alldftwocond <- alldf[, c(seq_len(9), idxcoltwoconds)]

        ## Reset column names so that .buildcolnames() in tepr() assigns them
        ## correctly for the two-condition subset
        colnames(alldftwocond) <- paste0("V", seq_len(ncol(alldftwocond)))

        ## Calling tepr on the defined conditions
        restepr <- .restepr(saveobjectpath, compname, reload, expdftwocond,
            alldftwocond, expthres, nbcpu, rounding, condonename, condtwoname,
            replaceval, pval, significant, windsizethres, countnathres,
            meancondonethres, meancondtwothres, pvaltheorythres,
            auccondonethreshigher, auccondonethreslower, auccondtwothres,
            attenuatedpvalksthres, outgrouppvalksthres, showtime, verbose)

        rm(alldftwocond)
        if (showmemory) print(gc()) else invisible(gc())
        return(restepr)

    }, verbose, expdf, alldf, expthres, nbcpu, rounding, replaceval, pval,
        significant, windsizethres, countnathres, meancondonethres,
        meancondtwothres, pvaltheorythres, auccondonethreshigher,
        auccondonethreslower, auccondtwothres, attenuatedpvalksthres,
        outgrouppvalksthres, saveobjectpath, reload, showtime, showmemory,
        simplify = FALSE)

    ## Naming each element with the comparison title
    names(reslist) <- apply(matcond, 2, function(currentcol) {
        return(paste(currentcol[1], currentcol[2], sep = "_vs_"))})

    if (!is.na(saveobjectpath)) {
        if (verbose) message("\t\t Saving to ", filepathname)
        saveRDS(reslist, file = filepathname)
    }

    return(reslist)
}


#' Perform tepr differential nascent rna-seq analysis for multiple conditions
#'
#' @description
#' This function performs tepr differential nascent RNA-seq analysis for
#' multiple conditions. It iterates over all pairwise comparisons of
#' conditions within the experiment table and calls the `tepr` function for
#' each pair.
#'
#' @usage
#' teprmulti(expdf, alldf, expthres, nbcpu = 1, rounding = 10,
#'  dontcompare = NULL, replaceval = NA, pval = 0.1, significant = FALSE,
#'  windsizethres = 50, countnathres = 20, pvaltheorythres = 0.1,
#'  meancondonethres = 0.5, meancondtwothres = 0.5,
#'  auccondonethreshigher = -10, auccondonethreslower = 15, auccondtwothres = 15,
#'  attenuatedpvalksthres = 2, outgrouppvalksthres = 0.2,
#'  saveobjectpath = NA, reload = FALSE, showtime = FALSE, showmemory = FALSE,
#'  verbose = TRUE)
#'
#' @param expdf A data frame containing experiment data that should have
#'          columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param alldf A data frame containing all transcript-related information,
#'          including biotype, chromosome, coordinates, transcript, gene,
#'          strand, window, ID and scores retrieved from the bedgraph files.
#' @param expthres A numeric value specifying the expression threshold.
#'          Transcripts with average expression values below this threshold
#'          will be filtered out from the returned transcript vector.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'          parallel computation on transcripts. Defaults to \code{1}.
#' @param rounding An integer specifying the rounding factor for computing ECDF.
#'          Default is \code{10}.
#' @param dontcompare An optional vector specifying conditions to exclude
#'          from the comparison. It should use condition names from expdf and
#'          follow the pattern \code{cond1_vs_cond2}. Defaults to \code{NULL}.
#' @param replaceval A value to replace non-significant attenuation values.
#'          Defaults to \code{NA}.
#' @param pval A numeric value specifying the p-value threshold for significance
#'          of the KS test. Defaults to \code{0.1}.
#' @param significant A logical indicating whether to filter out non-significant
#'          attenuation values. Defaults to \code{FALSE}.
#' @param windsizethres A numeric threshold for the minimum window size.
#'          Default is 50.
#' @param countnathres A numeric threshold for the maximum number of missing
#'          data points for an experiment (NA values). Default is 20.
#' @param pvaltheorythres A numeric threshold for the minimum p-value used to
#'          define the universe of genes. Default is 0.1.
#' @param meancondonethres A numeric threshold for the minimum mean transcription
#'          value in the first condition of each comparison. Default is 0.5.
#' @param meancondtwothres A numeric threshold for the minimum mean transcription
#'          value in the second condition of each comparison. Default is 0.5.
#' @param auccondonethreshigher A numeric threshold for the lower bound of the
#'          first condition AUC value in the outgroup classification.
#'          Default is -10.
#' @param auccondonethreslower A numeric threshold for the upper bound of the
#'          first condition AUC value in the outgroup classification.
#'          Default is 15.
#' @param auccondtwothres A numeric threshold for the minimum second condition
#'          AUC value used to classify attenuated genes. Default is 15.
#' @param attenuatedpvalksthres A numeric threshold for the negative log10 of
#'          the p-value (from KS test) for defining attenuated genes.
#'          Default is 2.
#' @param outgrouppvalksthres A numeric threshold for the maximum KS p-value
#'          used to define the outgroup. Default is 0.2.
#' @param saveobjectpath A character string specifying the path to save
#'          the object of the results for each comparison. The file names are
#'          defined with the 'condition' column of expdf. Defaults to \code{NA}.
#' @param reload Logical. If `TRUE`, reloads existing saved objects to avoid
#'  recomputation. Default is `FALSE`. If the function failed during object
#'  saving, make sure to delete the corresponding object.
#' @param showtime A logical value indicating if the duration of the function
#'          processing should be indicated before ending. Defaults to
#'          \code{FALSE}.
#' @param showmemory A logical value indicating if memory usage should be
#'          printed during the execution. Defaults to \code{FALSE}.
#' @param verbose A logical flag indicating whether to print progress messages.
#'          Defaults to \code{TRUE}.
#'
#' @return
#' A list of lists. Each inner list corresponds to a pairwise comparison
#' and contains two data frames:
#'   - `resmeandiff_<comparison>`: Results from the `meandifference` function
#'     for the specific comparison.
#'   - `resunigroupatt_<comparison>`: Results from the `universegroup` function
#'     for the specific comparison.
#'
#' @examples
#' ## Supposing the data have more than one condition
#' \dontrun{
#'   exptabpath <- "exp.csv"
#'   alldfpath <- "result-preprocessing.tsv"
#'   expdf <- read.csv(exptabpath)
#'   alldf <- read.delim(alldfpath, header = FALSE)
#'   expthres <- 0.1
#'   reslist <- teprmulti(expdf, alldf, expthres)}
#'
#' @seealso
#' [tepr]
#'
#' @export

teprmulti <- function(expdf, alldf, expthres, nbcpu = 1, rounding = 10,
    dontcompare = NULL, replaceval = NA, pval = 0.1, significant = FALSE,
    windsizethres = 50, countnathres = 20, pvaltheorythres = 0.1,
    meancondonethres = 0.5, meancondtwothres = 0.5, auccondonethreshigher = -10,
    auccondonethreslower = 15, auccondtwothres = 15, attenuatedpvalksthres = 2,
    outgrouppvalksthres = 0.2, saveobjectpath = NA, reload = FALSE,
    showtime = FALSE, showmemory = FALSE, verbose = TRUE) {

    if (showtime) start_teprmulti <- Sys.time()

    if (!length(unique(expdf$condition)) > 2)
        stop("\n[tepr] Error: Too few conditions.\n",
            "  Found 2 or fewer conditions. Use tepr() instead.\n")

    checkexptab(expdf)

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Building col names for alldf
    alldf <- .buildcolnames(expdf, alldf)

    ## Eliminating comparisons if dontcompare not NULL
    matcond <- .dontcompare(dontcompare, expdf, verbose)

    ## Calling tepr by pairs of contions
    filepathname <- file.path(saveobjectpath, "allcomplist.rds")

    if (!reload || !file.exists(filepathname)) {
        reslist <- .reslist(filepathname, matcond, verbose, expdf, alldf,
            expthres, nbcpu, rounding, replaceval, pval, significant,
            windsizethres, countnathres, meancondonethres, meancondtwothres,
            pvaltheorythres, auccondonethreshigher, auccondonethreslower,
            auccondtwothres, attenuatedpvalksthres, outgrouppvalksthres,
            saveobjectpath, reload, showtime, showmemory)
    } else {
        if (verbose) message("\t\t\t Loading ", filepathname)
        reslist <- readRDS(filepathname)
    }

    if (showtime) {
      end_teprmulti <- Sys.time()
      timing <- end_teprmulti - start_teprmulti
      message("\t\t ## Analysis teprmulti performed in: ",
        format(timing, digits = 2))
    }

    return(reslist)
}
