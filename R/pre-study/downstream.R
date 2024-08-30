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


##################
# PARAMETERS
##################


alldfpath <- "/g/romebioinfo/tmp/preprocessing/completeframedf.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
expthres <- 0.1
## Parallelization on bedgraph files. The maximum should be equal to the number of bedgraph files.  # nolint
nbcpubg <- 8
## Parallelization on transcripts. The maximum should be limited to the capacity of your machine.  # nolint
nbcputrans <- 20


##################
#FUNCTIONS
##################

averageandfilterexprs <- function(expdf, alldf, expthres, verbose = FALSE) { # nolint

    scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction,
      "score")

    ## Calculate the average expression per transcript (over each frame)
    if (verbose) message("\t Calculating average expression per transcript") # nolint
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint

    ## Remove a line if it contains only values < expthres (separating strands)
    if (verbose) message("\t Removing lines with values < expthres") # nolint
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
            stop("The element ", xname, # nolint
                " should be unique, contact the developer.") # nolint
}

.createecdfmat <- function(scoremat, rounding, transtable, direction) {

  ecdfmat <- apply(scoremat, 2, function(x, rounding, coordvec) {
    extendedcoordvec <- rep(coordvec, ceiling(x * rounding))
    fx <- ecdf(extendedcoordvec)(coordvec)
    return(fx)
  }, rounding, transtable$coord, simplify = TRUE)
  colnames(ecdfmat) <- gsub(direction, "", colnames(ecdfmat))
  colnames(ecdfmat) <- paste("Fx", colnames(ecdfmat), sep = "_") # nolint
  return(ecdfmat)
}

.checkdirection <- function(str, transtable) {

  if (isTRUE(all.equal(str, "-"))) {
    transtable <- transtable[order(transtable$coord), ]
    directionfill <- "updown"
  } else {
    directionfill <- "downup"
  }
  return(list(transtable, directionfill))
}

.computeecdf <- function(transtable, expdf, rounding, colnamevec) { # nolint

        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        .checkunique(str, "str")
        colnamestr <- colnamevec[which(expdf$strand == str)]

        ## If the strand is negative, re-order by coordinates
        direclist <- .checkdirection(str, transtable)
        transtable <- direclist[[1]]
        directionfill <- direclist[[2]]

        ## Building a matrix containing only the scores in the right direction
        ## for each experiment. Filling the NA values with tidyr::fill.
        scoremat <- transtable[, colnamestr]
        scoremat <- scoremat %>% tidyr::fill(contains("score"), # nolint
          .direction = directionfill)

        ## Replace the scores of transtable with the filled one
        transtable[, colnamestr] <- scoremat

        ## Retrieving the direction (fwd or rev) according to the transcript
        ## strand. It will be used to change the column names of scoremat and
        ## remove columns from transtab.
        direction <- unique(expdf[which(expdf$strand == str), "direction"])
        .checkunique(direction, "direction")
        opposedirect <- unique(expdf[which(expdf$strand != str), "direction"])
        .checkunique(opposedirect, "opposite direction") # nolint

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- .createecdfmat(scoremat, rounding, transtable, direction)

        ## Remove opposite strand from transtable and erase strand substring
        transtable <- transtable[, -grep(opposedirect, colnames(transtable))]
        colnames(transtable) <- gsub(direction, "", colnames(transtable))

        res <- cbind(transtable, ecdfmat)
        return(res)
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1, # nolint
  verbose = FALSE) {

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]

    ## Filtering the main table to keep only the expressed transcripts
    if (verbose) message("\t Filtering to keep only the expressed transcripts") # nolint
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    if (isTRUE(all.equal(length(idxnoexpr), 0)))
      warning("All the transcripts are expressed", immediate. = TRUE) # nolint
    else
      maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    if (verbose) message("\t Splitting the table by each transcript") # nolint
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")
    colnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction,
      "score")

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        colnamevec, rounding) {
        res <- .computeecdf(transtable, expdf, rounding, colnamevec)
        return(res)
    }, expdf, colnamevec, rounding, mc.cores = nbcpu)

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
          ## The function rowDiffs of the package matrixStats substracts the second argument to the first one. To respect the code just above, # nolint
          ## The indexes must be inverted with rev: meancolnames[rev(idxvec)]] -> for instance, given the two columns "mean_value_HS" and "mean_value_ctrl" # nolint
          ## as input, the function rowDiffs will do the subtraction "mean_value_ctrl" - "mean_value_HS" # nolint
          res <- matrixStats::rowDiffs(as.matrix(
            resmean[, meancolnames[rev(idxvec)]]))
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

.meandiffscorefx <- function(idxcondlist, df, nbwindows, currentcond,
    colnamevec, verbose) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, nbwindows,
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
                diffres <- meandf - df$coord / nbwindows
                colnames(diffres) <- paste0("diff_", idxname, "_", currentcond)
                res <- cbind(meandf, diffres)
            } else {
                res <- meandf
            }
            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, nbwindows,
            currentcond, colnamevec, verbose), SIMPLIFY = FALSE)

        return(meandifflist)
}

createmeandiff <- function(resultsecdf, expdf, nbwindows, verbose = FALSE) {

    ## for each condition, creates three columns:
    ##   - "mean_value_ctrl", "mean_Fx_ctrl", "diff_Fx_ctrl"
    ##   - "mean_value_HS", "mean_Fx_HS", "diff_Fx_HS"
    condvec <- unique(expdf$condition)
    rescondlist <- lapply(condvec, function(currentcond, df, nbwindows,
      verbose) {

        if (verbose) message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)
        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxscorefx(df, idxcond)

        ## The difference is used to calculate the AUC later on
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, nbwindows,
            currentcond, colnamevec, verbose)
        names(meandifflist) <- NULL
        meandiffres <- do.call("cbind", meandifflist)
        return(meandiffres)
    }, resultsecdf, nbwindows, verbose)
    resmean <- do.call("cbind", rescondlist)
    ## Computing all differences on mean columns
    if (verbose) message("Commputing all differences on mean columns")
    matdiff <- .creatematdiff(condvec, resmean)

    res <- cbind(resmean, matdiff)
    if (!isTRUE(all.equal(nrow(resultsecdf), nrow(res))))
        stop("The results of mean and diff should have the same number of ",
            "rows than resultsecdf, contact the developer")

    return(cbind(resultsecdf, res))
}


.returninfodf <- function(transtab, nbwindows = NULL) {

  transcript <- unique(transtab$transcript)
  gene <- unique(transtab$gene)
  strand <- unique(transtab$strand)
        .checkunique(transcript, "transcript-dauc_allconditions")
        .checkunique(gene, "gene-dauc_allconditions")
        .checkunique(strand, "strand-dauc_allconditions")
        infodf <- data.frame(transcript, gene, strand)

        if (!is.null(nbwindows)) {
            if (isTRUE(all.equal(as.character(strand), "+"))) {
                windsize <- floor(
                    (transtab$end[nbwindows] - transtab$start[1]) / nbwindows)
            } else {
                windsize <- floor(
                    (transtab$end[1] - transtab$start[nbwindows]) / nbwindows)
            }
            infodf <- cbind(infodf, windsize)
        }
        return(infodf)
}

dauc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1,
    dontcompare = NULL) {

    condvec <- unique(expdf$condition)
    resdflist <- mclapply(bytranslist, function(transtab, condvec) {

        ## Retrieve the column names for each comparison
        idxctrl <- grep("ctrl", condvec) # Cannot be empty, see checkexptab
        name1 <- paste0("mean_Fx_", condvec[idxctrl])
        name2 <- paste0("mean_Fx_", condvec[-idxctrl])
        diffname <- paste0("Diff_meanFx_", condvec[-idxctrl], "_",
          condvec[idxctrl])

        ## Perform a kolmogorov-smirnoff test between the two columns
        resks <- suppressWarnings(ks.test(transtab[, name1], transtab[, name2]))

        ## Calculate the area under the curve of the difference of means
        ## -> delta AUC
        deltadauc <- pracma::trapz(transtab[, "coord"], transtab[, diffname])
        ## Retrieve the p-value
        pvaldeltadaucks <- resks$p.value
        ## The KS test statistic is defined as the maximum value of the
        ## difference between A and B’s cumulative distribution functions (CDF)
        statdeltadaucks <- resks$statistic

        ## Build a one line data.frame with the proper col names
        ksdeltadaucdf <- data.frame(deltadauc, pvaldeltadaucks, statdeltadaucks)
        colnames(ksdeltadaucdf) <- paste(colnames(ksdeltadaucdf), name2,
            sep = "_")

        ## Retrieving transcript information
        infodf <- .returninfodf(transtab, nbwindows)

        ## Combining the two df as result
        resdf <- cbind(infodf, ksdeltadaucdf)
        return(resdf)
    }, condvec, mc.cores = nbcpu)

    resdf <- do.call("rbind", resdflist)
    return(resdf)
}



.buildaucdf <- function(transtab, difffxname, resks, meanvalname,
  currentcond) {
    auc <- pracma::trapz(transtab[, "coord"], transtab[, difffxname])
    pvalaucks <- resks$p.value
    stataucks <- resks$statistic
    meanvaluefull <- mean(transtab[, meanvalname])
    aucdf <- data.frame(auc, pvalaucks, stataucks, meanvaluefull)
    colnames(aucdf) <- paste(colnames(aucdf), currentcond, sep = "_")
    rownames(aucdf) <- paste(.returninfodf(transtab), collapse = "-")
    transinfo <- data.frame(transcript = transtab[1, "transcript"],
                    gene = transtab[1, "gene"], strand = transtab[1, "strand"])
    aucdf <- cbind(transinfo, aucdf)
    return(aucdf)
}

auc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1) {

  cumulative <- seq(1, nbwindows) / nbwindows
  condvec <- unique(expdf$condition)

  resdflist <- mclapply(bytranslist, function(transtab, condvec, cumulative,
    nbwindows) {
      ## Computing AUC, pval, and stat for each condition
      resauclist <- lapply(condvec, function(currentcond, transtab,
        cumulative) {
          ## Definition of column names
          difffxname <- paste0("diff_Fx_", currentcond)
          meanvalname <- paste0("mean_value_", currentcond)
          meanfxname <- paste0("mean_Fx_", currentcond)

          ## Perform a kolmogorov-smirnoff test between mean_Fx and cum.density
          resks <- suppressWarnings(ks.test(transtab[, meanfxname], cumulative))
          ## Build data.frame with auc information for the current transcript
          aucdf <- .buildaucdf(transtab, difffxname, resks, meanvalname,
            currentcond)
          return(aucdf)
      }, transtab, cumulative)
      aucdf <- do.call("cbind", resauclist)
      return(aucdf)
  }, condvec, cumulative, nbwindows, mc.cores = nbcpu)

  aucallconditions <- do.call("rbind", resdflist)
  return(aucallconditions)
}



countna <- function(allexprsdfs, expdf, nbcpu, verbose = FALSE) {

  maintable <- allexprsdfs[[1]]
  scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction,
    "score")
  condvec <- unique(expdf$condition)

  ## Splitting the table by each transcript
  if (verbose) message("\t Splitting the table by each transcript") # nolint
  transdflist <- split(maintable, factor(maintable$transcript))

  ## For each transcript
  nabytranslist <- parallel::mclapply(transdflist,
    function(transtable, scorecolvec, condvec) {
        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        colnamestr <- scorecolvec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]

        ## Counting NA for each condition (c=condition, m=matrix, n=colnames)
        res <- lapply(condvec, function(c, m, n) {
          length(which(apply(m[, grep(c, n)], 1, function(x) all(is.na(x)))))
        }, scoremat, colnamestr)
        res <- do.call("cbind", res)
        colnames(res) <- paste0(condvec, "_NA")
        info <- data.frame(gene = unique(transtable$gene),
          transcript = unique(transtable$transcript), strand = str)
        return(cbind(info, res))
    }, scorecolvec, condvec, mc.cores = nbcpu)
  return(do.call("rbind", nabytranslist))
}



.retrievekneeandmax <- function(condvec, transtable) { # nolint

  reslist <- lapply(condvec, function(cond, transtable) {

    difffxname <- paste0("diff_Fx_", cond)
    difffxvec <- transtable[, difffxname]

    ## If equality of difference within the same gene it takes the closest
    ## knee from the TSS # nolint
    resrow <- transtable[which(difffxvec == max(difffxvec)), ] %>% # nolint
          dplyr::slice_min(coord, n = 1) # nolint
    res <- data.frame(resrow$coord, resrow[, difffxname])
    colnames(res) <- c(paste0("knee_AUC_", cond), paste0("max_", difffxname))

    return(res)
      }, transtable)
    return(reslist)
}

kneeid <- function(transdflist, expdf, nbcputrans, verbose = FALSE) {

  condvec <- unique(expdf$condition)

  bytransres <- parallel::mclapply(transdflist, function(transtable, condvec) {
      bycondreslist <- .retrievekneeandmax(condvec, transtable)
       return(cbind(transcript = transtable$transcript[1],
        do.call("cbind", bycondreslist)))
    }, condvec, mc.cores = nbcputrans)
  res <- do.call("rbind", bytransres)
  return(res)
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
message("Filtering transcripts based on expression") # nolint
allexprsdfs <- averageandfilterexprs(expdf, alldf, expthres)
message("Calculating ECDF") # nolint
start_time <- Sys.time()
resecdf <- genesECDF(allexprsdfs, expdf, nbcpu = nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
resultsecdf <- resecdf[[1]]
nbwindows <- resecdf[[2]]
saveRDS(resultsecdf, "/g/romebioinfo/tmp/downstream/resultsecdf.rds")

start_time <- Sys.time()
message("Calculating means and differences")
dfmeandiff <- createmeandiff(resultsecdf, expdf, nbwindows)
saveRDS(dfmeandiff, "/g/romebioinfo/tmp/downstream/dfmeandiff.rds")
## Splitting result by transcripts
message("\t Splitting results by transcripts")
bytranslistmean <- split(dfmeandiff, factor(dfmeandiff$transcript))
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint

message("Computing the differences (d or delta) of AUC")
start_time <- Sys.time()
dfaucallcond <- dauc_allconditions(bytranslistmean, expdf, nbwindows,
  nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
saveRDS(dfaucallcond, "/g/romebioinfo/tmp/downstream/dfaucallcond.rds")

# Calculate the Area Under Curve (AUC), All conditions vs y=x
# Calculate Mean Value over the full gene body in All conditions.
message("Computing the Area Under Curve (AUC)")
start_time <- Sys.time()
aucallcond <- auc_allconditions(bytranslistmean, expdf, nbwindows,
  nbcpu = nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
saveRDS(aucallcond, "/g/romebioinfo/tmp/downstream/aucallcond.rds")

message("Calculating number of missing values for each transcript and for",
  " each condition")
start_time <- Sys.time()
matnatrans <- countna(allexprsdfs, expdf, nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
saveRDS(matnatrans, "/g/romebioinfo/tmp/downstream/matnatrans.rds")

message("Retrieving knee and max")
start_time <- Sys.time()
kneedf <- kneeid(bytranslistmean, expdf, nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
saveRDS(kneedf, "/g/romebioinfo/tmp/downstream/kneedf.rds")



!!!!!!!!!!!!!!!!!
Attenuation_fun <- function(AUC_KS_Knee_NA_DF, concat_df, pval,Replaced) {

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

  Complete_summary <-  left_join(concat_df, AUC_KS_Knee_NA_DF,
    by = c("gene","transcript","strand"))

  for (cond in Conditions) {
    mean_value_condi_name <- paste0("mean_value_", cond)
    print(mean_value_condi_name)
    knee_column_name <- paste0("knee_AUC_", cond)
    Attenuation_cond <- paste0("Attenuation_", cond)
    UPmean_cond <- paste0("UP_mean_", cond)
    DOWNmean_cond <- paste0("DOWN_mean_", cond)
    AUC_KS_Knee_NA_DF[[Attenuation_cond]] <- NA

    result <- Complete_summary %>% group_by(transcript) %>% arrange(coord) %>%
      dplyr::reframe(transcript=transcript[1],
        !!sym(UPmean_cond) := mean((!!sym(mean_value_condi_name))[coord <= !!sym(knee_column_name)]), # nolint
        !!sym(DOWNmean_cond) := mean((!!sym(mean_value_condi_name))[coord >= !!sym(knee_column_name) & coord <= max(coord)])) %>% # nolint
        select(transcript,!!sym(UPmean_cond),!!sym(DOWNmean_cond), !!sym(DOWNmean_cond)) %>% distinct() # nolint

    AUC_KS_Knee_NA_DF <- left_join(AUC_KS_Knee_NA_DF,result, by=c("transcript"))
    AUC_KS_Knee_NA_DF[[Attenuation_cond]] <- 100 - AUC_KS_Knee_NA_DF[[DOWNmean_cond]]/AUC_KS_Knee_NA_DF[[UPmean_cond]]*100 # nolint
  }

  if (exists("Replaced") && !is.na(Replaced)) {
    if (Replaced != "NOT") {
      for (cond in Conditions) {
        p_AUC_cond <- paste0("p_AUC_", cond)
        print(p_AUC_cond)
        AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
          mutate(!!paste0("Attenuation_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, # nolint
          Replaced, .data[[paste0("Attenuation_", cond)]])) ## replacing the Attenuation by an inout value is KS test > at threshold # nolint
        AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
          mutate(!!paste0("knee_AUC_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, .data[[paste0("knee_AUC_", cond)]])) ## replacing the knee by NA is KS test > at threshold # nolint
      }
    }
  } else {
    for (cond in Conditions) {
      p_AUC_cond <- paste0("p_AUC_", cond)
      print(p_AUC_cond)
      AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
        mutate(!!paste0("Attenuation_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, # nolint
        .data[[paste0("Attenuation_", cond)]])) ## replacing the Attenuation by an input value if KS test > at threshold # nolint
      AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
        mutate(!!paste0("knee_AUC_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, # nolint
        .data[[paste0("knee_AUC_", cond)]])) ## replacing the knee by NA if KS test > at threshold # nolint
    }
  }
  return(AUC_KS_Knee_NA_DF)
}

!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!! SUMMARY IN ONE TABLE OF ALL THE VALUES COMPUTED ABOVE
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

!!!!!!!!!!!!!!! THIS ENABLES A FILTERING ON NA, WINDOWSIZE, ETC
!!!!!!!!!!!!!!!!!! SEE IF CAN BE INTEGRATED SOMEWHERE
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