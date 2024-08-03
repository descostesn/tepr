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
library("matrixStats")



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

averageandfilterexprs <- function(expdf, alldf, expthres, verbose = FALSE) { # nolint

    scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ## Calculate the average expression per transcript (over each frame)
    if(verbose) message("\t Calculating average expression per transcript") # nolint
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint
    ## Remove a line if it contains only values < expthres (separating strands)
    if(verbose) message("\t Removing lines with values < expthres")
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

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1,
  verbose = FALSE) { # nolint

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]

    ## Filtering the main table to keep only the expressed transcripts
    if (verbose) message("\t Filtering to keep only the expressed transcripts")
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    if (verbose) message("\t Splitting the table by each transcript")
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")
    framevec <- seq_len(nbrows)
    colnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
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
    colnamevec, verbose) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, tosub, nbrows,
            currentcond, colnamevec, verbose) {
            if (verbose)
              message("\t Calculating average and difference between ",
                "replicates for columns '", idxname, "' of ", currentcond)

            ## Calculating the column of mean scores for currentcond
            ## The result is a data.frame made of a single column
            if (length(idxvalvec) >= 2)
                meandf <- data.frame(rowMeans(df[, idxvalvec], na.rm = FALSE))
            else
                meandf <- df[, idxvalvec]
            colnames(meandf) <- paste0("mean_", idxname, "_", currentcond)

            if (isTRUE(all.equal(idxname, "Fx"))) {
                diffres <- meandf - tosub
                colnames(diffres) <- paste0("diff_", idxname, "_", currentcond)
                res <- cbind(meandf, diffres)
            } else {
                res <- meandf
            }
            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, tosub, nbrows,
            currentcond, colnamevec, verbose), SIMPLIFY = FALSE)

        return(meandifflist)
}

createmeandiff <- function(resultsecdf, expdf, verbose = FALSE) {

    ## for each condition, creates three columns:
    ##   - "mean_value_ctrl", "mean_Fx_ctrl", "diff_Fx_ctrl"
    ##   - "mean_value_HS", "mean_Fx_HS", "diff_Fx_HS"
    condvec <- unique(expdf$condition)
    rescondlist <- lapply(condvec, function(currentcond, df,
      verbose) {

        if (verbose) message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)

        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxscorefx(df, idxcond)

        ## The difference is used to calculate the AUC later on
        nbrows <- nrow(df)
        tosub <- df$window / nbrows
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, tosub, nbrows,
            currentcond, colnamevec, verbose)
        names(meandifflist) <- NULL

        meandiffres <- do.call("cbind", meandifflist)

        return(meandiffres)
    }, resultsecdf, verbose)

    resmean <- do.call("cbind", rescondlist)

    ## Computing all differences on mean columns
    ##   - "Diff_meanValue_ctrl_HS", "Diff_meanValue_HS_ctrl"
    ##   - "Diff_meanFx_ctrl_HS", "Diff_meanFx_HS_ctrl"
    if (verbose) message("Commputing all differences on mean columns")
    .creatematdifflist
    categoryvec <- c("value", "Fx")
    matdifflist <- lapply(categoryvec, function(currentcat, condvec, resmean) {
      meancolnames <- paste("mean", currentcat, condvec, sep = "_")
      ## Generating all combinations of elements (combn not good)
      idxvec <- seq_len(length(condvec))
      matidx <- matrix(c(idxvec, rev(idxvec)), ncol = 2)
      difflist <- apply(matidx, 2, function(idxvec, meancolnames, resmean,
        currentcat, condvec) {
        res <- matrixStats::rowDiffs(as.matrix(resmean[,meancolnames[idxvec]]))
        colnamestr <- paste("Diff", paste0("mean", currentcat),
          paste(condvec[idxvec], collapse = "_"), sep = "_")
        res <- as.vector(res)
        attr(res, "name") <- colnamestr
        return(res)
      }, meancolnames, resmean, currentcat, condvec, simplify = FALSE)
      diffmat <- do.call("cbind", difflist)
      colnames(diffmat) <- sapply(difflist, function(x) attributes(x)$name)
      return(diffmat)
    }, condvec, resmean)
    matdiff <- do.call("cbind", matdifflist)

!!
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
message("Filtering transcripts based on expression")
allexprsdfs <- averageandfilterexprs(expdf, alldf, expthres)
message("Calculating ECDF")
resultsecdf <- genesECDF(allexprsdfs, expdf, nbcpu = nbcpu)
dfmeandiff <- createmeandiff(resultsecdf, expdf)










> head(concat_dfFX_res,2)
         biotype  chr     coor1     coor2         transcript gene strand window
1 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +      1
2 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +      2
                           id         ctrl_rep1         ctrl_rep2
1 ENST00000000233.10_ARF5_+_1 ctrl_rep1.forward ctrl_rep2.forward
2 ENST00000000233.10_ARF5_+_2 ctrl_rep1.forward ctrl_rep2.forward
          HS_rep1         HS_rep2 coord value_ctrl_rep1_score
1 HS_rep1.forward HS_rep2.forward     1                     0
2 HS_rep1.forward HS_rep2.forward     2                     0
  value_ctrl_rep2_score value_HS_rep1_score value_HS_rep2_score
1                     0                   0                   0
2                     0                   0                   0
  Fx_ctrl_rep1_score Fx_ctrl_rep2_score Fx_HS_rep1_score Fx_HS_rep2_score
1                  0                  0                0                0
2                  0                  0                0                0
  mean_value_ctrl mean_Fx_ctrl diff_Fx_ctrl mean_value_HS mean_Fx_HS diff_Fx_HS
1               0            0       -0.005             0          0     -0.005
2               0            0       -0.010             0          0     -0.010



> head(concat_Diff_mean_res,2)
         biotype  chr     coor1     coor2         transcript gene strand window
1 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +      1
2 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +      2
                           id         ctrl_rep1         ctrl_rep2
1 ENST00000000233.10_ARF5_+_1 ctrl_rep1.forward ctrl_rep2.forward
2 ENST00000000233.10_ARF5_+_2 ctrl_rep1.forward ctrl_rep2.forward
          HS_rep1         HS_rep2 coord value_ctrl_rep1_score
1 HS_rep1.forward HS_rep2.forward     1                     0
2 HS_rep1.forward HS_rep2.forward     2                     0
  value_ctrl_rep2_score value_HS_rep1_score value_HS_rep2_score
1                     0                   0                   0
2                     0                   0                   0
  Fx_ctrl_rep1_score Fx_ctrl_rep2_score Fx_HS_rep1_score Fx_HS_rep2_score
1                  0                  0                0                0
2                  0                  0                0                0
  mean_value_ctrl mean_Fx_ctrl diff_Fx_ctrl mean_value_HS mean_Fx_HS diff_Fx_HS
1               0            0       -0.005             0          0     -0.005
2               0            0       -0.010             0          0     -0.010


> head(dAUC_allcondi_res,2)
          transcript gene strand window_size
1 ENST00000000233.10 ARF5      +          16
2  ENST00000000412.8 M6PR      -          46



          transcript    gene strand window_size
1 ENST00000000233.10    ARF5      +          16
2  ENST00000000412.8    M6PR      -          46
3 ENST00000000442.11   ESRRA      +          56
4  ENST00000001008.6   FKBP4      +          52
5  ENST00000001146.7 CYP26B1      -          93
6  ENST00000002125.9 NDUFAF7      +          87




> head(count_NA_res)
# A tibble: 6 × 4
  gene    transcript         strand Count_NA
  <chr>   <chr>              <chr>     <int>
1 ARF5    ENST00000000233.10 +             0
2 ESRRA   ENST00000000442.11 +             0
3 FKBP4   ENST00000001008.6  +             0
4 NDUFAF7 ENST00000002125.9  +             0
5 SEMA3F  ENST00000002829.8  +             0
6 CFTR    ENST00000003084.11 +             2


> head(KneeID_res)
          transcript
1 ENST00000000233.10
2  ENST00000000412.8
3 ENST00000000442.11
4  ENST00000001008.6
5  ENST00000001146.7
6  ENST00000002125.9

> head(AUC_KS_Knee_NA.df)
          transcript    gene strand window_size Count_NA
1 ENST00000000233.10    ARF5      +          16        0
2  ENST00000000412.8    M6PR      -          46       16
3 ENST00000000442.11   ESRRA      +          56        0
4  ENST00000001008.6   FKBP4      +          52        0
5  ENST00000001146.7 CYP26B1      -          93        0
6  ENST00000002125.9 NDUFAF7      +          87        0


> head(AUC_KS_Knee_NA.df)
# A tibble: 6 × 9
  transcript         chr    coor1  coor2 strand gene   size window_size Count_NA
  <chr>              <chr>  <int>  <int> <chr>  <chr> <dbl>       <int>    <int>
1 ENST00000000233.10 chr7  1.28e8 1.28e8 +      ARF5   3290          16        0
2 ENST00000000412.8  chr12 8.94e6 8.95e6 -      M6PR   9285          46       16
3 ENST00000000442.11 chr11 6.43e7 6.43e7 +      ESRRA 11220          56        0
4 ENST00000001008.6  chr12 2.79e6 2.81e6 +      FKBP4 10454          52        0
5 ENST00000001146.7  chr2  7.21e7 7.21e7 -      CYP2… 18625          93        0
6 ENST00000002125.9  chr2  3.72e7 3.72e7 +      NDUF… 17503          87        0


> head(tst_df)
# A tibble: 6 × 9
  transcript         chr    coor1  coor2 strand gene   size window_size Count_NA
  <chr>              <chr>  <int>  <int> <chr>  <chr> <dbl>       <int>    <int>
1 ENST00000000233.10 chr7  1.28e8 1.28e8 +      ARF5   3290          16        0
2 ENST00000000412.8  chr12 8.94e6 8.95e6 -      M6PR   9285          46       16
3 ENST00000000442.11 chr11 6.43e7 6.43e7 +      ESRRA 11220          56        0
4 ENST00000001008.6  chr12 2.79e6 2.81e6 +      FKBP4 10454          52        0
5 ENST00000001146.7  chr2  7.21e7 7.21e7 -      CYP2… 18625          93        0
6 ENST00000002125.9  chr2  3.72e7 3.72e7 +      NDUF… 17503          87        0


> tst_df <- tst_df %>%
  mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 &
    !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
    !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
  relocate(Universe, .before = 1)
Error in `mutate()`:
ℹ In argument: `Universe = ifelse(...)`.
Caused by error:
! object 'MeanValueFull_ctrl' not found
Run `rlang::last_trace()` to see where the error occurred.


> tst_df <- tst_df %>% mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest)) >1.5, "Attenuated", NA), # nolint
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group) # nolint
  ) %>% relocate(Group, .before = 2)
Error in `mutate()`:
ℹ In argument: `Group = ifelse(...)`.
Caused by error:
! object 'Universe' not found
Run `rlang::last_trace()` to see where the error occurred.

