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

    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        framevec, colnamevec) {
        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        .checkunique(str, "str")
        colnamestr <- colnamevec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- apply(scoremat, 2, function(x, rounding, framevec) {
            extendedframevec <- rep(framevec, ceiling(x * rounding))
            fx <- ecdf(extendedframevec)(framevec)
            return(fx)
        }, rounding, framevec, simplify = TRUE)
        colnames(ecdfmat) <- paste("Fx", colnames(ecdfmat), sep = "_")
        res <- cbind(transtable, ecdfmat)
        return(res)
    }, expdf, framevec, colnamevec, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(concatdf)
}

.condcolidx <- function(currentcond, df) {
    idxcond <- grep(currentcond, colnames(df))
    if (isTRUE(all.equal(length(idxcond), 0)))
        stop("Problem in function meananddiff, condition not found in ",
                "column names. Contact the developer.")
    return(idxcond)
}

.idxcondlist <- function(df, idxcond) {
    idxcondfx <- grep("Fx", colnames(df[idxcond]))
    if (isTRUE(all.equal(length(idxcondfx), 0)))
        stop("Problem in function meananddiff, column Fx not found in ",
                "column names. Contact the developer.")
    idxcondlist <- list(scores = idxcond[-idxcondfx],
            fx = idxcond[idxcondfx])
    return(idxcondlist)
}

.meandiffscorefx <- function(idxcondlist, df, tosub, nbrows, currentcond,
    colnamevec) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, tosub, nbrows,
            condname, colnamevec){
            message("\t Calculating average and difference between replicates",
                " for columns ", idxname, " of ", condname)

            if (length(idxvalvec) >= 2)
                meanvec <- rowMeans(df[, idxvalvec], na.rm = FALSE)
            else
                meanvec <- df[, idxvalvec]

            diffres <- data.frame(df[, idxvalvec] - tosub) ## data.frame is used in case there is only one column # nolint
            colnames(diffres) <- paste(colnamevec[idxvalvec], "diff", sep = "_")
            meanname <- paste0("mean_", idxname, "_", condname)
            res <- cbind(meanvec, diffres)
            colnames(res)[1] <- meanname
            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, tosub, nbrows,
            currentcond, colnamevec), SIMPLIFY = FALSE)
        return(meandifflist)
}

meananddiff <- function(resultsecdf, expdf) {

    rescondlist <- lapply(unique(expdf$condition), function(currentcond, df) {

        message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)

        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxcondlist(df, idxcond)

        ## The difference is used to calculate the AUC
        nbrows <- nrow(df)
        tosub <- df$window / nbrows
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, tosub, nbrows,
            currentcond, colnamevec)

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
dfmeandiff <- meananddiff(resultsecdf, expdf)








> head(resultsECDF)
         biotype  chr     coor1     coor2         transcript gene strand window
1 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +      1
2 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +      2
3 protein-coding chr7 127588443 127588459 ENST00000000233.10 ARF5      +      3
4 protein-coding chr7 127588459 127588475 ENST00000000233.10 ARF5      +      4
5 protein-coding chr7 127588475 127588491 ENST00000000233.10 ARF5      +      5
6 protein-coding chr7 127588491 127588507 ENST00000000233.10 ARF5      +      6
                           id         ctrl_rep1         ctrl_rep2
1 ENST00000000233.10_ARF5_+_1 ctrl_rep1.forward ctrl_rep2.forward
2 ENST00000000233.10_ARF5_+_2 ctrl_rep1.forward ctrl_rep2.forward
3 ENST00000000233.10_ARF5_+_3 ctrl_rep1.forward ctrl_rep2.forward
4 ENST00000000233.10_ARF5_+_4 ctrl_rep1.forward ctrl_rep2.forward
5 ENST00000000233.10_ARF5_+_5 ctrl_rep1.forward ctrl_rep2.forward
6 ENST00000000233.10_ARF5_+_6 ctrl_rep1.forward ctrl_rep2.forward
          HS_rep1         HS_rep2 coord value_ctrl_rep1_score
1 HS_rep1.forward HS_rep2.forward     1              0.000000
2 HS_rep1.forward HS_rep2.forward     2              0.000000
3 HS_rep1.forward HS_rep2.forward     3              0.000000
4 HS_rep1.forward HS_rep2.forward     4              0.000000
5 HS_rep1.forward HS_rep2.forward     5              0.000000
6 HS_rep1.forward HS_rep2.forward     6              0.440169
  value_ctrl_rep2_score value_HS_rep1_score value_HS_rep2_score
1                     0                   0                   0
2                     0                   0                   0
3                     0                   0                   0
4                     0                   0                   0
5                     0                   0                   0
6                     0                   0                   0
  Fx_ctrl_rep1_score Fx_ctrl_rep2_score Fx_HS_rep1_score Fx_HS_rep2_score
1        0.000000000                  0                0                0
2        0.000000000                  0                0                0
3        0.000000000                  0                0                0
4        0.000000000                  0                0                0
5        0.000000000                  0                0                0
6        0.000497265                  0                0                0





> head(concat_dfFX_res)
         biotype  chr     coor1     coor2         transcript gene strand window
1 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +      1
2 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +      2
3 protein-coding chr7 127588443 127588459 ENST00000000233.10 ARF5      +      3
4 protein-coding chr7 127588459 127588475 ENST00000000233.10 ARF5      +      4
5 protein-coding chr7 127588475 127588491 ENST00000000233.10 ARF5      +      5
6 protein-coding chr7 127588491 127588507 ENST00000000233.10 ARF5      +      6
                           id         ctrl_rep1         ctrl_rep2
1 ENST00000000233.10_ARF5_+_1 ctrl_rep1.forward ctrl_rep2.forward
2 ENST00000000233.10_ARF5_+_2 ctrl_rep1.forward ctrl_rep2.forward
3 ENST00000000233.10_ARF5_+_3 ctrl_rep1.forward ctrl_rep2.forward
4 ENST00000000233.10_ARF5_+_4 ctrl_rep1.forward ctrl_rep2.forward
5 ENST00000000233.10_ARF5_+_5 ctrl_rep1.forward ctrl_rep2.forward
6 ENST00000000233.10_ARF5_+_6 ctrl_rep1.forward ctrl_rep2.forward
          HS_rep1         HS_rep2 coord value_ctrl_rep1_score
1 HS_rep1.forward HS_rep2.forward     1              0.000000
2 HS_rep1.forward HS_rep2.forward     2              0.000000
3 HS_rep1.forward HS_rep2.forward     3              0.000000
4 HS_rep1.forward HS_rep2.forward     4              0.000000
5 HS_rep1.forward HS_rep2.forward     5              0.000000
6 HS_rep1.forward HS_rep2.forward     6              0.440169
  value_ctrl_rep2_score value_HS_rep1_score value_HS_rep2_score
1                     0                   0                   0
2                     0                   0                   0
3                     0                   0                   0
4                     0                   0                   0
5                     0                   0                   0
6                     0                   0                   0
  Fx_ctrl_rep1_score Fx_ctrl_rep2_score Fx_HS_rep1_score Fx_HS_rep2_score
1        0.000000000                  0                0                0
2        0.000000000                  0                0                0
3        0.000000000                  0                0                0
4        0.000000000                  0                0                0
5        0.000000000                  0                0                0
6        0.000497265                  0                0                0