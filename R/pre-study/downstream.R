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
library("pracma")

!!!!!!!!!
environment: /g/romebioinfo/tmp/downstream
!!!!!!!!!!!


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

        ## Defining coordinates according to the strand
        if (isTRUE(all.equal(str, '+')))
            transtable <- cbind(transtable, coord = transtable$window)
        else
            transtable <- cbind(transtable, coord = rev(transtable$window))

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

    return(list(concatdf, nbrows))
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

.meandiffscorefx <- function(idxcondlist, df, nbrows, currentcond,
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
                diffres <- meandf - df$coord / nbrows
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

.creatematdiff <- function(condvec, resmean) {

  categoryvec <- c("value", "Fx")
  matdifflist <- lapply(categoryvec, function(currentcat, condvec, resmean) {
    meancolnames <- paste("mean", currentcat, condvec, sep = "_")

    ## Generating all combinations of elements (combn not good)
    idxvec <- seq_len(length(condvec))
    matidx <- matrix(c(idxvec, rev(idxvec)), ncol = 2)

    ## Generating differences of columns
    difflist <- apply(matidx, 2, function(idxvec, meancolnames, resmean,
        currentcat, condvec) {
          res <- matrixStats::rowDiffs(as.matrix(
            resmean[,meancolnames[idxvec]]))
          colnamestr <- paste("Diff", paste0("mean", currentcat),
            paste(condvec[idxvec], collapse = "_"), sep = "_")
          res <- as.vector(res)
          attr(res, "name") <- colnamestr
          return(res)
      }, meancolnames, resmean, currentcat, condvec, simplify = FALSE)

      ## Combining vectors into a matrix and defining col names
      diffmat <- do.call("cbind", difflist)
      colnames(diffmat) <- sapply(difflist, function(x) attributes(x)$name)
      return(diffmat)
  }, condvec, resmean)

  ## Building a matrix from the diff on values and Fx
  matdiff <- do.call("cbind", matdifflist)
  return(matdiff)
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
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, nbrows,
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
    matdiff <- .creatematdiff(condvec, resmean)

    ## Combining the matrices with mean columns and differences of the mean
    ## columns
    res <- cbind(resmean, matdiff)

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
resecdf <- genesECDF(allexprsdfs, expdf, nbcpu = nbcpu)
resultsecdf <- resecdf[[1]]
nbwindows <- resecdf[[2]]

message("Calculating means and differences")
dfmeandiff <- createmeandiff(resultsecdf, expdf)


!!!!!!!!!!!!!!!!!!!!!!!!!

# df=dfmeandiff
# dontcompare = NULL
dAUC_allcondi_fun <- function(df, expdf, nbwindows, dontcompare = NULL) {

  bytranslist <- split(df, factor(df$transcript))
  condvec <- unique(expdf$condition)
  resdflist <- lapply(bytranslist, function(transtab, condvec) {

    ## Sorting table according to strand
    transtab <- transtab[order(as.numeric(transtab$coord)), ]

    ## Retrieve the column names for each comparison
    idxctrl <- grep("ctrl", condvec) # Cannot be empty, see checkexptab
    name1 <- paste0("mean_Fx_", condvec[idxctrl])
    name2 <- paste0("mean_Fx_", condvec[-idxctrl])
    diffname <- paste0("Diff_meanFx_", condvec[-idxctrl], "_", condvec[idxctrl])

    ## Perform a kolmogorov-smirnoff test between the two columns
    resks <- suppressWarnings(ks.test(transtab[, name1], transtab[, name2]))

    ## Calculate the area under the curve of the difference of means - delta AUC
    deltaauc <- pracma::trapz(transtab[,"coord"], transtab[, diffname])
    ## Retrieve the p-value
    pvalks <- resks$p.value
    ## The KS test statistic is defined as the maximum value of the difference
    ## between A and B’s cumulative distribution functions (CDF)
    statks <- resks$statistic

    ## Build a one line data.frame with the proper col names
    ksaucdf <- data.frame(deltaauc, pvalks, statks)
    colnames(ksaucdf) <- paste(colnames(ksaucdf), name2, sep = "_")

    ## Retrieving transcript information
    transcript <- unique(transtab$transcript)
    gene <- unique(transtab$gene)
    strand <- unique(transtab$strand)
    .checkunique(transcript, "transcript-dAUC_allcondi_fun")
    .checkunique(gene, "gene-dAUC_allcondi_fun")
    .checkunique(strand, "strand-dAUC_allcondi_fun")
    if (isTRUE(all.equal(strand, '+')))
        windsize <- floor(
            (transtab$end[nbwindows] - transtab$start[1])/nbwindows)
    else
        windsize <- floor(
            (transtab$end[1] - transtab$start[nbwindows])/nbwindows)
    infodf <- data.frame(transcript, gene, strand, windsize)

    ## Combining the two df as result
    resdf <- cbind(infodf, ksaucdf)
    return(resdf)
  }, condvec)

!! CURRENT !!  resdf <- do.call("rbind", resdflist)

  ## Create a data.frame with the columns: transcript, gene, strand, window_size
  !!!!!!!!!!!!!!!!!!! This pre-creation could be avoided with a split
  dAUC_allcondi <- df  %>% 
    filter(window==round(nbwindows/2))  %>%
    mutate(window_size = abs(coor2-coor1), .keep = "all") %>%
    select("transcript", "gene", "strand", "window_size") %>% distinct()

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

!!!!!!!!!!!!!!!!!!!!!! checkexptab in preprocessing.R ensure a maximum of two condition
!!!!!!!!!!!!!!!!!!!!!!!!!! The for loops are not necessary
  for (i in 1:length(Conditions)) {
    for (j in (1+i):length(Conditions)) {
      if (j > length(Conditions)){ break}
      cond1 <- Conditions[i]
      cond2 <- Conditions[j]

      ## making sure to not do useless comparison for the user,
      compare <- paste0(cond1, " vs ", cond2)
      if (! compare %in% dontcompare) {
        mean_Fx_condi_name1 <- paste0("mean_Fx_", cond1) #e,g Control
        mean_Fx_condi_name2 <- paste0("mean_Fx_", cond2) #e.g HS

        Diff_meanFx_name1 <- paste0("Diff_meanFx_",cond1,"_",cond2)
        Diff_meanFx_name2 <- paste0("Diff_meanFx_",cond2,"_",cond1) 

        df_name <- paste0("gene_summary_AUC_", Diff_meanFx_name1)
        assign(df_name,data.frame()) 

        # Get the data frame using the dynamically generated name
        !!!!!!!!!!!!!
        >  dAUC
[1] "dAUC_Diff_meanFx_HS_ctrl"
>  p_dAUC
[1] "p_dAUC_Diff_meanFx_HS_ctrl"
> D_dAUC
[1] "D_dAUC_Diff_meanFx_HS_ctrl"
!!!!!!!!!!!!!!!
        dAUC_summary_condi_df <- data.frame()
        dAUC_summary_condi_df <- get(df_name)
        dAUC <- paste0("dAUC_",Diff_meanFx_name2)
        p_dAUC <- paste0("p_dAUC_",Diff_meanFx_name2)
        D_dAUC <- paste0("D_dAUC_",Diff_meanFx_name2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! This performs a treatment by transcript, use a split
        dAUC_summary_condi_df  <- df %>%
          group_by(transcript) %>% arrange(coord) %>%
          reframe(gene=gene[1], strand=strand, transcript=transcript, 
            !!dAUC := trapz(coord,!!sym(Diff_meanFx_name2)),# delta AUC
            !!p_dAUC := ks.test(!!sym(mean_Fx_condi_name1),
            !!sym(mean_Fx_condi_name2))$p.value,
            !!D_dAUC := ks.test(!!sym(mean_Fx_condi_name1),
            !!sym(mean_Fx_condi_name2))$statistic,) %>% dplyr::distinct()
        
        !!!!!!!!!!!!!
        > head( as.data.frame(dAUC_summary_condi_df))
          transcript    gene strand dAUC_Diff_meanFx_HS_ctrl
1 ENST00000000233.10    ARF5      +                7.3183914
2  ENST00000000412.8    M6PR      -               -0.6347988
3 ENST00000000442.11   ESRRA      +                4.5515080
4  ENST00000001008.6   FKBP4      +                2.7953480
5  ENST00000001146.7 CYP26B1      -               -6.1229193
6  ENST00000002125.9 NDUFAF7      +                0.3563332
  p_dAUC_Diff_meanFx_HS_ctrl D_dAUC_Diff_meanFx_HS_ctrl
1                  0.1122497                      0.120
2                  1.0000000                      0.020
3                  0.8642828                      0.060
4                  0.4653198                      0.085
5                  0.2202056                      0.105
6                  1.0000000                      0.025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        assign(df_name, dAUC_summary_condi_df)
        dAUC_allcondi <- left_join(dAUC_allcondi, dAUC_summary_condi_df, by = c("transcript", "gene","strand")) # nolint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
> head(dAUC_allcondi)
          transcript    gene strand window_size
1 ENST00000000233.10    ARF5      +          16
2  ENST00000000412.8    M6PR      -          46
3 ENST00000000442.11   ESRRA      +          56
4  ENST00000001008.6   FKBP4      +          52
5  ENST00000001146.7 CYP26B1      -          93
6  ENST00000002125.9 NDUFAF7      +          87
> dAUC_allcondi <- left_join(dAUC_allcondi, dAUC_summary_condi_df, by = c("transcript", "gene","strand")) # nolint
> head(dAUC_allcondi)
          transcript    gene strand window_size dAUC_Diff_meanFx_HS_ctrl
1 ENST00000000233.10    ARF5      +          16                7.3183914
2  ENST00000000412.8    M6PR      -          46               -0.6347988
3 ENST00000000442.11   ESRRA      +          56                4.5515080
4  ENST00000001008.6   FKBP4      +          52                2.7953480
5  ENST00000001146.7 CYP26B1      -          93               -6.1229193
6  ENST00000002125.9 NDUFAF7      +          87                0.3563332
  p_dAUC_Diff_meanFx_HS_ctrl D_dAUC_Diff_meanFx_HS_ctrl
1                  0.1122497                      0.120
2                  1.0000000                      0.020
3                  0.8642828                      0.060
4                  0.4653198                      0.085
5                  0.2202056                      0.105
6                  1.0000000                      0.025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      }
    }
  }

  dAUC_allcondi <- dAUC_allcondi %>%
  mutate(across(contains("p_dAUC"), ~ modify_p_values(.))) 

  # Get the column names containing "p_dAUC"
  p_dAUC_columns <- grep("p_dAUC", colnames(dAUC_allcondi), value = TRUE)

  # Loop through each column, calculate adjusted p-values, and create new columns
  for (col_name in p_dAUC_columns) {
    print(col_name)
    adjusted_col_name <- paste0("adjFDR_", col_name)
    dAUC_allcondi[[adjusted_col_name]] <- p.adjust(dAUC_allcondi[[col_name]], method = "fdr")
  }
  return(dAUC_allcondi)
}


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


> head(concat_Diff_mean_res)
         biotype  chr     coor1     coor2         transcript gene strand window
1 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +      1
2 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +      2
                           id         ctrl_rep1         ctrl_rep2
1 ENST00000000233.10_ARF5_+_1 ctrl_rep1.forward ctrl_rep2.forward
2 ENST00000000233.10_ARF5_+_2 ctrl_rep1.forward ctrl_rep2.forward
          HS_rep1         HS_rep2 coord value_ctrl_rep1_score
1 HS_rep1.forward HS_rep2.forward     1              0.000000
2 HS_rep1.forward HS_rep2.forward     2              0.000000
  value_ctrl_rep2_score value_HS_rep1_score value_HS_rep2_score
1                     0                   0                   0
2                     0                   0                   0
  Fx_ctrl_rep1_score Fx_ctrl_rep2_score Fx_HS_rep1_score Fx_HS_rep2_score
1        0.000000000                  0                0                0
2        0.000000000                  0                0                0
  mean_value_ctrl mean_Fx_ctrl diff_Fx_ctrl mean_value_HS mean_Fx_HS diff_Fx_HS
1       0.0000000 0.0000000000  -0.00500000             0          0     -0.005
2       0.0000000 0.0000000000  -0.01000000             0          0     -0.010
  Diff_meanValue_ctrl_HS Diff_meanValue_HS_ctrl Diff_meanFx_ctrl_HS
1              0.0000000              0.0000000        0.0000000000
2              0.0000000              0.0000000        0.0000000000
  Diff_meanFx_HS_ctrl
1        0.0000000000
2        0.0000000000


head(dfmeandiff,2)
         biotype  chr     start       end         transcript gene strand window
1 protein_coding chr7 127588411 127588426 ENST00000000233.10 ARF5      +      1
2 protein_coding chr7 127588427 127588442 ENST00000000233.10 ARF5      +      2
                           id ctrl1_score ctrl2_score HS1_score HS2_score
1 ENST00000000233.10_ARF5_+_1           0           0         0         0
2 ENST00000000233.10_ARF5_+_2           0           0         0         0
  Fx_ctrl1_score Fx_ctrl2_score Fx_HS1_score Fx_HS2_score mean_value_ctrl
1              0              0            0            0               0
2              0              0            0            0               0
  mean_Fx_ctrl  diff_Fx_ctrl mean_value_HS mean_Fx_HS    diff_Fx_HS
1            0 -3.310381e-07             0          0 -3.310381e-07
2            0 -6.620763e-07             0          0 -6.620763e-07
  Diff_meanvalue_ctrl_HS Diff_meanvalue_HS_ctrl Diff_meanFx_ctrl_HS
1                      0                      0                   0
2                      0                      0                   0
  Diff_meanFx_HS_ctrl
1                   0
2                   0





Time difference of 30.59541 secs
> head(dAUC_allcondi_res,2)
          transcript gene strand window_size dAUC_Diff_meanFx_HS_ctrl
1 ENST00000000233.10 ARF5      +          16                7.3183914
2  ENST00000000412.8 M6PR      -          46               -0.6347988
  p_dAUC_Diff_meanFx_HS_ctrl D_dAUC_Diff_meanFx_HS_ctrl
1                  0.1122497                       0.12
2                  1.0000000                       0.02
  adjFDR_p_dAUC_Diff_meanFx_HS_ctrl
1                         0.2732393
2                         1.0000000


Time difference of 1.000898 mins
> head(AUC_allcondi_res,2)
          transcript gene strand window_size    AUC_ctrl p_AUC_ctrl D_AUC_ctrl
1 ENST00000000233.10 ARF5      +          16 -16.1578100 0.01195204      0.160
2  ENST00000000412.8 M6PR      -          46   0.4324419 0.99999997      0.025
  MeanValueFull_ctrl    AUC_HS  p_AUC_HS D_AUC_HS MeanValueFull_HS
1           4.877027 -8.839419 0.3274975    0.095         4.564402
2           9.180490 -0.202357 0.9996971    0.035         9.789388
  adjFDR_p_AUC_ctrl adjFDR_p_AUC_HS
1        0.04014826       0.4433998
2        1.00000000       1.0000000

Time difference of 5.873777 secs
> head(count_NA_res,2)
# A tibble: 2 × 4
  gene  transcript         strand Count_NA
  <chr> <chr>              <chr>     <int>
1 ARF5  ENST00000000233.10 +             0
2 ESRRA ENST00000000442.11 +             0


Time difference of 6.39923 secs
> head(KneeID_res,2)
          transcript knee_AUC_ctrl max_diff_Fx_ctrl knee_AUC_HS max_diff_Fx_HS
1 ENST00000000233.10           200       0.00000000         115     0.01909246
2  ENST00000000412.8            29       0.02321981          29     0.02434586




> AUC_KS_Knee_NA.df <- left_join(AUC_allcondi_res, dAUC_allcondi_res,
  by = c("transcript", "gene", "strand", "window_size"))  %>%
  left_join(., KneeID_res, by = c("transcript"))  %>%
  left_join(., count_NA_res, by = c("gene", "transcript", "strand"))
> AUC_KS_Knee_NA.df <- concat_Diff_mean_res %>% group_by(transcript) %>%
  summarise( chr=chr[1], coor1=min(coor1), coor2=max(coor2), strand=strand[1],
  gene=gene[1], transcript=transcript[1], size=coor2-coor1+1) %>%
  left_join(AUC_KS_Knee_NA.df, by=c("gene", "transcript", "strand"))
> head(AUC_KS_Knee_NA.df,2)
# A tibble: 2 × 27
  transcript         chr    coor1  coor2 strand gene   size window_size AUC_ctrl
  <chr>              <chr>  <int>  <int> <chr>  <chr> <dbl>       <int>    <dbl>
1 ENST00000000233.10 chr7  1.28e8 1.28e8 +      ARF5   3290          16  -16.2
2 ENST00000000412.8  chr12 8.94e6 8.95e6 -      M6PR   9285          46    0.432
# ℹ 18 more variables: p_AUC_ctrl <dbl>, D_AUC_ctrl <dbl>,
#   MeanValueFull_ctrl <dbl>, AUC_HS <dbl>, p_AUC_HS <dbl>, D_AUC_HS <dbl>,
#   MeanValueFull_HS <dbl>, adjFDR_p_AUC_ctrl <dbl>, adjFDR_p_AUC_HS <dbl>,
#   dAUC_Diff_meanFx_HS_ctrl <dbl>, p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   D_dAUC_Diff_meanFx_HS_ctrl <dbl>, adjFDR_p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   knee_AUC_ctrl <dbl>, max_diff_Fx_ctrl <dbl>, knee_AUC_HS <dbl>,
#   max_diff_Fx_HS <dbl>, Count_NA <int>


> head(tst_df,2)
# A tibble: 2 × 33
  transcript         chr    coor1  coor2 strand gene   size window_size AUC_ctrl
  <chr>              <chr>  <int>  <int> <chr>  <chr> <dbl>       <int>    <dbl>
1 ENST00000000233.10 chr7  1.28e8 1.28e8 +      ARF5   3290          16  -16.2
2 ENST00000000412.8  chr12 8.94e6 8.95e6 -      M6PR   9285          46    0.432
# ℹ 24 more variables: p_AUC_ctrl <dbl>, D_AUC_ctrl <dbl>,
#   MeanValueFull_ctrl <dbl>, AUC_HS <dbl>, p_AUC_HS <dbl>, D_AUC_HS <dbl>,
#   MeanValueFull_HS <dbl>, adjFDR_p_AUC_ctrl <dbl>, adjFDR_p_AUC_HS <dbl>,
#   dAUC_Diff_meanFx_HS_ctrl <dbl>, p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   D_dAUC_Diff_meanFx_HS_ctrl <dbl>, adjFDR_p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   knee_AUC_ctrl <dbl>, max_diff_Fx_ctrl <dbl>, knee_AUC_HS <dbl>,
#   max_diff_Fx_HS <dbl>, Count_NA <int>, Attenuation_ctrl <dbl>, …




> mean_value_control_full <- "MeanValueFull_ctrl"
mean_value_stress <- "MeanValueFull_HS"
AUC_ctrl <- "AUC_ctrl"
AUC_stress <- "AUC_HS"
p_value_KStest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
p_value_theoritical<- "adjFDR_p_AUC_ctrl"
tst_df <- tst_df %>%
  mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 &
    !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
    !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
  relocate(Universe, .before = 1)
tst_df <- tst_df %>% mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest)) >1.5, "Attenuated", NA), # nolint
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group) # nolint
  ) %>% relocate(Group, .before = 2)
> head(tst_df)
# A tibble: 6 × 35
  Universe Group   transcript chr    coor1  coor2 strand gene   size window_size
  <lgl>    <chr>   <chr>      <chr>  <int>  <int> <chr>  <chr> <dbl>       <int>
1 FALSE    NA      ENST00000… chr7  1.28e8 1.28e8 +      ARF5   3290          16
2 FALSE    NA      ENST00000… chr12 8.94e6 8.95e6 -      M6PR   9285          46
3 TRUE     Outgro… ENST00000… chr11 6.43e7 6.43e7 +      ESRRA 11220          56
4 TRUE     Outgro… ENST00000… chr12 2.79e6 2.81e6 +      FKBP4 10454          52
5 FALSE    NA      ENST00000… chr2  7.21e7 7.21e7 -      CYP2… 18625          93
6 TRUE     Outgro… ENST00000… chr2  3.72e7 3.72e7 +      NDUF… 17503          87
# ℹ 25 more variables: AUC_ctrl <dbl>, p_AUC_ctrl <dbl>, D_AUC_ctrl <dbl>,
#   MeanValueFull_ctrl <dbl>, AUC_HS <dbl>, p_AUC_HS <dbl>, D_AUC_HS <dbl>,
#   MeanValueFull_HS <dbl>, adjFDR_p_AUC_ctrl <dbl>, adjFDR_p_AUC_HS <dbl>,
#   dAUC_Diff_meanFx_HS_ctrl <dbl>, p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   D_dAUC_Diff_meanFx_HS_ctrl <dbl>, adjFDR_p_dAUC_Diff_meanFx_HS_ctrl <dbl>,
#   knee_AUC_ctrl <dbl>, max_diff_Fx_ctrl <dbl>, knee_AUC_HS <dbl>,
#   max_diff_Fx_HS <dbl>, Count_NA <int>, Attenuation_ctrl <dbl>, …