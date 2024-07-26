####################################
# This script goes through documentation/explore.R and homogenizes it with
# preprocessing.R
#
# Descostes - R-4.4.1 - July 2024
####################################

library("tidyr")
library("dplyr")
library("tidyselect")
library("parallel")


##################
# PARAMETERS
##################

alldfpath <- "/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv"
expthres <- 0.1
nbcpu <- 10


##################
#FUNCTIONS
##################

averageandfilterexprs <- function(expdf, alldf, expthres) { # nolint

    scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ## Calculate the average expression per transcript (over each frame)
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint
    ## Remove a line if it contains only values < expthres (separating strands)
    dfstrandlist <- mapply(function(strandname, directname, dfbytrans,
        expthres) {
            if ((isTRUE(all.equal(strandname, "+")) &&
                isTRUE(all.equal(directname, "rev"))) ||
                (isTRUE(all.equal(strandname, "-")) &&
                isTRUE(all.equal(directname, "fwd"))))
                    stop("Strand and direction do not match, contact the ",
                        "developer")
            dfstrand <- dfbytranscript %>%
                dplyr::filter(strand == strandname) %>% # nolint
                dplyr::select(gene, transcript, strand, # nolint
                tidyselect::contains(directname))  %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ !is.na(.))) %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ . > expthres))
            return(dfstrand)
        }, unique(expdf$strand), unique(expdf$direction),
            MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)

    exptranstab <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>%
        dplyr::arrange(transcript) %>% dplyr::pull(transcript) # nolint

    return(list(maintable = alldf, exptranlist = exptranstab))
}



.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname,
                " should be unique, contact the developer.")
}

.computeecdf <- function(transtable, expdf, rounding, framevec, colnamevec) {
        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        .checkunique(str, "str")
        colnamestr <- colnamevec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]
        direction <- unique(expdf[which(expdf$strand == str), "direction"])
        .checkunique(direction, "direction")
        opposedirect <- unique(expdf[which(expdf$strand != str), "direction"])
        .checkunique(opposedirect, "opposite direction")

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- apply(scoremat, 2, function(x, rounding, framevec) {
            extendedframevec <- rep(framevec, ceiling(x * rounding))
            fx <- ecdf(extendedframevec)(framevec)
            return(fx)
        }, rounding, framevec, simplify = TRUE)
        colnames(ecdfmat) <- gsub(direction, "", colnames(ecdfmat))
        colnames(ecdfmat) <- paste("Fx", colnames(ecdfmat), "score", sep = "_")

        ## Remove opposite strand from transtable and erase strand substring
        transtable <- transtable[,-grep(opposedirect, colnames(transtable))]
        colnames(transtable) <- gsub(direction, "_score", colnames(transtable))

        res <- cbind(transtable, ecdfmat)
        return(res)
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1) { # nolint

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]

    ## Filtering the main table to keep only the expressed transcripts
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")
    framevec <- seq_len(nbrows)
    colnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ## Computing ecdf on each transcript
    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        framevec, colnamevec, rounding) {
        res <- .computeecdf(transtable, expdf, rounding, framevec, colnamevec)
        return(res)
    }, expdf, framevec, colnamevec, rounding, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(concatdf)
}



.condcolidx <- function(currentcond, df) {
    idxcond <- grep(currentcond, colnames(df))
    if (isTRUE(all.equal(length(idxcond), 0)))
        stop("Problem in function createmeandiff, condition not found in ",
                "column names. Contact the developer.")
    return(idxcond)
}

.idxscorefx <- function(df, idxcond) {
    idxcondfx <- grep("Fx", colnames(df[idxcond]))
    if (isTRUE(all.equal(length(idxcondfx), 0)))
        stop("Problem in function createmeandiff, column Fx not found in ",
                "column names. Contact the developer.")
    idxcondlist <- list(value = idxcond[-idxcondfx],
            Fx = idxcond[idxcondfx])
    return(idxcondlist)
}

.meandiffscorefx <- function(idxcondlist, df, tosub, nbrows, currentcond,
    colnamevec) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, tosub, nbrows,
            condname, colnamevec){
            message("\t Calculating average and difference between replicates",
                " for columns '", idxname, "' of ", condname)

            ## Calculating the column of mean scores for currentcond
            ## The result is a data.frame made of a single column
            if (length(idxvalvec) >= 2)
                meandf <- data.frame(rowMeans(df[, idxvalvec], na.rm = FALSE))
            else
                meandf <- df[, idxvalvec]
            colnames(meandf) <- paste0("mean_", idxname, "_", condname)

            if (isTRUE(all.equal(idxname, "Fx"))) {
                diffres <- meandf - tosub
                colnames(diffres) <- paste0("diff_", idxname, "_", condname)
                res <- cbind(meandf, diffres)
            } else {
                res <- meandf
            }

            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, tosub, nbrows,
            currentcond, colnamevec), SIMPLIFY = FALSE)
        return(meandifflist)
}

createmeandiff <- function(resultsecdf, expdf) {

    ## for each condition, creates three columns:
    ##   - "mean_value_ctrl", "mean_Fx_ctrl", "diff_Fx_ctrl"
    ##   - "mean_value_HS", "mean_Fx_HS", "diff_Fx_HS"

    rescondlist <- lapply(unique(expdf$condition), function(currentcond, df) {

        message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)

        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxscorefx(df, idxcond)

        ## The difference is used to calculate the AUC later on
        nbrows <- nrow(df)
        tosub <- df$window / nbrows
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, tosub, nbrows,
            currentcond, colnamevec)
        names(meandifflist) <- NULL

        meandiffres <- do.call("cbind", meandifflist)
        return(meandiffres)
    }, resultsecdf)

    res <- do.call("cbind", rescondlist)
    if (!isTRUE(all.equal(nrow(resultsecdf), nrow(res))))
        stop("The results of mean and diff should have the same number of ",
            "rows than resultsecdf, contact the developer")

    return(cbind(resultsecdf, res))
}


##################
# MAIN
##################

## Reading alldf and info tab
alldf <- readRDS(alldfpath)
expdf <- read.csv(exptabpath, header = TRUE)

## Filtering out non expressed transcripts:
## 1) for each column, calculate the average expression per transcript (over each frame) # nolint
## 2) For each column, remove a line if it contains only values < expthres separating strands # nolint
allexprsdfs <- averageandfilterexprs(expdf, alldf, expthres)
resultsecdf <- genesECDF(allexprsdfs, expdf, nbcpu = nbcpu)
dfmeandiff <- createmeandiff(resultsecdf, expdf)
