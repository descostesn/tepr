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
filtertabpath <- "/g/romebioinfo/Projects/tepr/Dataset/filtertab.csv"

expthres <- 0.1
## Parallelization on bedgraph files. The maximum should be equal to the number of bedgraph files.  # nolint
nbcpubg <- 8
## Parallelization on transcripts. The maximum should be limited to the capacity of your machine.  # nolint
nbcputrans <- 20
testonerep <- FALSE


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
            if (verbose) {
              message("\t Calculating average and difference between ",
                "replicates for columns '", idxname, "' of ", currentcond)
              if (isTRUE(all.equal(length(idxvalvec), 1)))
                warning("Only one replicate, copy scores to mean columns",
                  immediate. = TRUE)
            }

            ## Calculating the column of mean scores for currentcond
            ## The result is a data.frame made of a single column
            if (length(idxvalvec) >= 2) {
                meandf <- data.frame(rowMeans(df[, idxvalvec], na.rm = FALSE))
            } else {
                meandf <- as.data.frame(df[, idxvalvec])
            }
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

.dauc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1,
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

    ## Correct p-values using FDR
    idx <- grep("pvaldeltadaucks", colnames(resdf))
    fdrvec <- p.adjust(resdf[, idx], method = "fdr")

    resdf <- cbind(resdf, fdrvec)
    colnamevec <- colnames(resdf)
    idxfdr <- which(colnamevec == "fdrvec")
    colnames(resdf)[idxfdr] <- paste0("adjFDR_", colnamevec[idx])
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

.auc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1) {

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
  idxdup <- which(duplicated(colnames(aucallconditions)))
  aucallconditions <- aucallconditions[, -idxdup]

  ## Correcting p-val with FDR
  idxpvalvec <- grep("pvalaucks", colnames(aucallconditions))
  fdrlist <- lapply(idxpvalvec, function(idxpval, tab) {
    return(p.adjust(tab[, idxpval], method = "fdr"))
  }, aucallconditions)
  fdrdf <- do.call("cbind", fdrlist)
  colnames(fdrdf) <- paste0("adjFDR_", colnames(aucallconditions)[idxpvalvec]) # nolint

  aucallconditions <- cbind(aucallconditions, fdrdf)
  return(aucallconditions)
}

allauc <- function(bytranslistmean, expdf, nbwindows, nbcputrans,
  dontcompare = NULL, verbose = TRUE) {

    if (verbose) message("\t Computing the differences (d or delta) of AUC")
    start_time <- Sys.time()
    daucallcond <- .dauc_allconditions(bytranslistmean, expdf, nbwindows,
      nbcputrans)
    end_time <- Sys.time()
    if (verbose) message("\t\t ## Analysis performed in: ",
      end_time - start_time) # nolint
    #saveRDS(dfaucallcond, "/g/romebioinfo/tmp/downstream/dfaucallcond.rds") # nolint

    ## Calculate the Area Under Curve (AUC), All conditions vs y=x
    ## Calculate Mean Value over the full gene body in All conditions.
    if (verbose) message("\t Computing the Area Under Curve (AUC)")
    start_time <- Sys.time()
    aucallcond <- .auc_allconditions(bytranslistmean, expdf, nbwindows,
      nbcpu = nbcputrans)
    end_time <- Sys.time()
    if (verbose) message("\t\t ## Analysis performed in: ",
      end_time - start_time) # nolint
    #saveRDS(aucallcond, "/g/romebioinfo/tmp/downstream/aucallcond.rds") # nolint

    ## Merging the two tables by transcript
    if (verbose) message("Merging results")
    allauc <- merge(aucallcond, daucallcond,
      by = c("gene", "transcript", "strand"))
    return(allauc)
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
        res <- sapply(condvec, function(c, m, n) {
          idxcol <- grep(c, n)
          if (!isTRUE(all.equal(length(idxcol), 1))) {
            return(length(which(apply(m[, idxcol], 1,
              function(x) all(is.na(x))))))
          } else { ## Only one replicate
            return(length(which(is.na(m[, idxcol]))))
          }
        }, scoremat, colnamestr)

        ## Retrieving total NA and transcript info
        countna <- sum(res)
        info <- data.frame(gene = unique(transtable$gene),
          transcript = unique(transtable$transcript), strand = str)
        return(cbind(info, countna))
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


.summarytrans <- function(bytransmeanlist, nbcpu) {
  summarydflist <- mclapply(bytranslistmean, function(trans) {
    coor1 <- min(trans$start)
    coor2 <- max(trans$end)
    return(data.frame(chr = trans$seqnames[1], coor1, coor2,
          strand = trans$strand[1], gene = trans$gene[1],
          transcript = trans$transcript[1], size = coor2 - coor1 + 1))
  }, mc.cores = nbcpu)
  summarydf <- do.call("rbind", summarydflist)
  return(summarydf)
}

.computeupdown <- function(completbytrans, condvec, nbcpu) {

  updownbytranslist <- mclapply(completbytrans, function(trans, condvec) {
    ## For each condition
    updownlist <- lapply(condvec, function(cond, trans) {
      kneecolname <- paste0("knee_AUC_", cond)
      meancolname <- paste0("mean_value_", cond)

      idxup <- which(trans$coord <= trans[, kneecolname])
      if (isTRUE(all.equal(length(idxup), 0)))
        stop("Problem in retrieving idxup, contact the developer.")
      upmean <- mean(trans[idxup, meancolname])

      idxdown <- which(trans$coord >= trans[, kneecolname] &
                          trans$coord <= max(trans$coord))
      if (isTRUE(all.equal(length(idxdown), 0)))
        stop("Problem in retrieving idxdown, contact the developer.")
      downmean <- mean(trans[idxdown, meancolname])

      ## Calculating attenuation
      att <- 100 - downmean / upmean * 100

      res <- data.frame(trans$transcript[1], upmean, downmean, att)
      colnames(res) <- c("transcript", paste0("upmean-", cond),
              paste0("downmean-", cond), paste0("attenuation-", cond))
      return(res)
    }, trans)

    return(do.call("cbind", updownlist))
  }, condvec, mc.cores = nbcpu)

  updowndf <- do.call("rbind", updownbytranslist)
  updowndf <- updowndf[, -which(duplicated(colnames(updowndf)))]
  return(updowndf)
}

attenuation <- function(allaucdf, kneedf, matnatrans, bytranslistmean, expdf,
  dfmeandiff, nbcpu = 1, verbose = TRUE) {

      if (verbose) message("\t Merging tables")
      allaucknee <- merge(allaucdf, kneedf, by = "transcript")
      mergecolnames <- c("gene", "transcript", "strand")
      allauckneena <- merge(allaucknee, matnatrans, by = mergecolnames)

      if (verbose) message("\t Building summary")
      summarydf <- .summarytrans(bytranslistmean, nbcpu)
      if (verbose) message("\t Merging summary")
      auckneenasum <- merge(summarydf, allauckneena, by = mergecolnames)

      ## Merging the mean table with the previous one
      if (verbose) message("\t Merging detailed mean table with summary")
      complet <- merge(dfmeandiff, auckneenasum, by = mergecolnames)

      ## Splitting the previous table by transcript
      if (verbose) message("\t Splitting the previous table by transcript")
      completbytrans <- split(complet, factor(complet$transcript))
      condvec <- unique(expdf$condition)

      ## For each transcript
      if (verbose) message("\t Computing up and down mean")
      updowndf <- .computeupdown(completbytrans, condvec, nbcpu)

      ## Merging attenuation to the complete table
      auckneenasumatt <- merge(auckneenasum, updowndf, by = "transcript")
      return(auckneenasumatt)
}


.createboolmat <- function(tab, completedf) {

  booleanlist <- apply(tab, 1, function(currentfilter, completedf) {

    currentfeature <- as.character(currentfilter["feature"])

    if (isTRUE(all.equal(currentfeature, "countna"))) {
      return(as.vector(completedf["countna"] < currentfilter["threshold"]))

    } else if (isTRUE(all.equal(currentfeature, "windowsize"))) {
      return(as.vector(completedf["windsize"] > currentfilter["threshold"]))

    } else if (isTRUE(all.equal(currentfeature, "fullmean"))) {
      colstr <- paste0("meanvaluefull_", currentfilter["condition"])
      return(completedf[, colstr] > currentfilter["threshold"])

    } else if (isTRUE(all.equal(currentfeature, "pvalauc"))) {
      colstr <- paste0("adjFDR_pvalaucks_", currentfilter["condition"]) # nolint
      return(completedf[, colstr] > currentfilter["threshold"])

    } else if (isTRUE(all.equal(currentfeature, "auc"))) {
      colstr <- paste0("auc_", currentfilter["condition"])
      if (!is.na(currentfilter["threshold2"]))
        res <- completedf[, colstr] > currentfilter["threshold"] &
          completedf[, colstr] < currentfilter["threshold2"]
      else
        res <- completedf[, colstr] > currentfilter["threshold"]
      return(res)
    } else if (isTRUE(all.equal(currentfeature, "daucfdrlog10"))) {
      colnamevec <- colnames(completedf)
      idx <- grep("adjFDR_pvaldeltadaucks_mean", colnamevec) # nolint
      .checkunique(idx)
      colstr <- colnamevec[idx]
      return(-log10(completedf[, colstr]) > currentfilter["threshold"])
    } else {
      stop("The feature ", currentfeature, " is not handlded.",
        "Allowed features are: countna, windowsize, fullmean, pvalauc, auc, ",
        "and daucfdrlog10. If you are sure that there is no typo, contact the",
        " developper.")
    }
  }, completedf, simplify = FALSE)

  booleanmat <- do.call("cbind", booleanlist)
  colnames(booleanmat) <- paste(tab[, "feature"], tab[, "condition"], sep = "-")
  return(booleanmat)
}

universegroup <- function(completedf, expdf, filterdf, verbose = TRUE) {

  ## Retrieving universe and group information
  universetab <- filterdf[which(filterdf$universe), ]
  grouptab <- filterdf[which(filterdf$group), ]

  ## Creating bool matrix and vector for universe
  universemat <- .createboolmat(universetab, completedf)
  universevec <- apply(universemat, 1, all)

  ## Creating bool matrix and vector for group
  groupvec <- rep(NA, nrow(completedf))
  groupmat <- .createboolmat(grouptab, completedf)
  attenuatedcols <-  as.vector(apply(filterdf[which(filterdf$group.attenuated),
    c("feature", "condition")], 1, paste, collapse = "-"))
  outgroupcols <-  as.vector(apply(filterdf[which(filterdf$group.outgroup),
    c("feature", "condition")], 1, paste, collapse = "-"))
  attenuatedvec <- apply(cbind(universevec, groupmat[, attenuatedcols]), 1, all)
  outgroupvec <- apply(cbind(universevec, groupmat[, outgroupcols]), 1, all)
  if (any(attenuatedvec))
    groupvec[attenuatedvec] <- "Attenuated"
  else
    warning("No attenuated genes were found")
  if (any(outgroupvec))
    groupvec[outgroupvec] <- "Outgroup"
  else
    warning("No outgroup genes were found")

  ## Building the final data.frame
  res <- data.frame(universe = universevec, group = groupvec)
  unigroupdf <- cbind(res, completedf)

  return(unigroupdf)
}


checkfilter <- function(filterdf, expdf) {

  filtercond <- unique(filterdf$condition[!is.na(filterdf$condition)])
  if (!isTRUE(all.equal(filtercond, unique(expdf$condition))))
    stop("The condition column of your experiment and filter tab should",
      " contain the same values.")

  if (all(!filterdf$universe))
    stop("All the rows of the universe column are set to FALSE. No rows will",
      " be used for the analysis.")

  if (all(!filterdf$group))
    stop("All the rows of the group column are set to FALSE. No rows will",
      " be used for the analysis.")

  if (all(!filterdf$group.attenuated))
    stop("All the rows of the group.attenuated column are set to FALSE. No ",
      "rows will be used for the analysis.")

    if (all(!filterdf$group.outgroup))
    stop("All the rows of the group.attenuated column are set to FALSE. No ",
      "rows will be used for the analysis.")

  featurevec <- c("auc", "countna", "windowsize", "fullmean", "daucfdrlog10",
    "pvalauc")
  idx <- match(featurevec, filterdf$feature)
  idxna <- which(is.na(idx))
  if (!isTRUE(all.equal(length(idxna), 0)))
    stop("The following features are missing from your filter table: ",
      paste(featurevec[idxna], collapse = "-"))

  colnametab <- colnames(filterdf)
  colnamevec <- c("condition", "feature", "threshold", "threshold2", "universe",
    "group", "group.attenuated", "group.outgroup")

  if (any(duplicated(colnametab)))
    stop("The columns of the filter table should be unique and must contain:",
      paste(colnamevec, collapse = "-"))

  idx <- match(colnamevec, colnametab)
  idxna <- which(is.na(idx))
  if (!isTRUE(all.equal(length(idxna), 0)))
    stop("The following columns are missing from your filter table: ",
      colnamevec[idxna])
}


##################
# MAIN
##################

## Reading alldf and info tab
alldf <- readRDS(alldfpath)
expdf <- read.csv(exptabpath, header = TRUE)
filterdf <- read.csv(filtertabpath, header = TRUE)
checkfilter(filterdf, expdf)

if (testonerep) {
  ## Removing the second replicate
  idxrep2 <- grep("ctrl2|HS2", colnames(alldf))
  alldf <- alldf[, -idxrep2]
  idxrep2 <- which(expdf$replicate == 2)
  expdf <- expdf[-idxrep2, ]
}

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
if (!testonerep) {
  saveRDS(resultsecdf, "/g/romebioinfo/tmp/downstream/resultsecdf.rds")
} else {
  saveRDS(resultsecdf, "/g/romebioinfo/tmp/downstream/resultsecdf-onerep.rds")
}

## Calculating means and differences
start_time <- Sys.time()
message("Calculating means and differences")
dfmeandiff <- createmeandiff(resultsecdf, expdf, nbwindows)
if (!testonerep) {
  saveRDS(dfmeandiff, "/g/romebioinfo/tmp/downstream/dfmeandiff.rds")
} else {
  saveRDS(dfmeandiff, "/g/romebioinfo/tmp/downstream/dfmeandiff-onerep.rds")
}
## Splitting result by transcripts
message("\t Splitting results by transcripts")
bytranslistmean <- split(dfmeandiff, factor(dfmeandiff$transcript))
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint

## Computing the differences (d or delta) of AUC and calculate the Area Under
## Curve (AUC), All conditions vs y=x
## Calculate Mean Value over the full gene body in All conditions.
message("AUC and differences")
allaucdf <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)
if (!testonerep) {
  saveRDS(allaucdf, "/g/romebioinfo/tmp/downstream/allaucdf.rds")
} else {
  saveRDS(allaucdf, "/g/romebioinfo/tmp/downstream/allaucdf-onerep.rds")
}

message("Calculating number of missing values for each transcript and for",
  " each condition")
start_time <- Sys.time()
matnatrans <- countna(allexprsdfs, expdf, nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
if (!testonerep) {
  saveRDS(matnatrans, "/g/romebioinfo/tmp/downstream/matnatrans.rds")
} else {
  saveRDS(matnatrans, "/g/romebioinfo/tmp/downstream/matnatrans-onerep.rds")
}


message("Retrieving knee and max")
start_time <- Sys.time()
kneedf <- kneeid(bytranslistmean, expdf, nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
if (!testonerep) {
  saveRDS(kneedf, "/g/romebioinfo/tmp/downstream/kneedf.rds")
} else {
  saveRDS(kneedf, "/g/romebioinfo/tmp/downstream/kneedf-onerep.rds")
}


message("Calculating attenuation values")
start_time <- Sys.time()
completedf <- attenuation(allaucdf, kneedf, matnatrans, bytranslistmean, expdf,
  dfmeandiff, nbcpu = nbcputrans)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
if (!testonerep) {
  saveRDS(completedf, "/g/romebioinfo/tmp/downstream/completedf.rds")
} else {
  saveRDS(completedf, "/g/romebioinfo/tmp/downstream/completedf-onerep.rds")
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
unigroupdf <- readRDS("/g/romebioinfo/tmp/downstream/unigroupdf.rds")
completedf <- readRDS("/g/romebioinfo/tmp/downstream/completedf.rds")
tst_df <- readRDS("/g/romebioinfo/tmp/explore/tst_df-1.rds")


####### With tst_df
mean_value_control_full <- "MeanValueFull_ctrl"
mean_value_stress <- "MeanValueFull_HS"
AUC_ctrl <- "AUC_ctrl"
AUC_stress <- "AUC_HS"
p_value_KStest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
p_value_theoritical<- "adjFDR_p_AUC_ctrl"

tstdf <- tst_df %>%
  mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 &
    !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
    !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
  relocate(Universe, .before = 1)

tstdf <- tstdf %>% mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest)) >1.5, "Attenuated", NA), # nolint
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group) # nolint
  ) %>% relocate(Group, .before = 2)

> print(table(tstdf$Universe))

FALSE  TRUE
 8373  6612
> print(table(tstdf$Group))

Attenuated   Outgroup
       513       5374

####### With completedf

mean_value_control_full <- "meanvaluefull_ctrl"
mean_value_stress <- "meanvaluefull_HS"
AUC_ctrl <- "auc_ctrl"
AUC_stress <- "auc_HS"
p_value_KStest <- "adjFDR_pvaldeltadaucks_mean_Fx_HS"
p_value_theoritical<- "adjFDR_pvalaucks_ctrl"

completetstdf <- completedf %>%
  mutate(Universe = ifelse(windsize > 50 & countna < 20 &
    !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
    !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
  relocate(Universe, .before = 1)

completetstdf <- completetstdf %>% mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest)) >1.5, "Attenuated", NA), # nolint
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group) # nolint
  ) %>% relocate(Group, .before = 2)

> print(table(completetstdf$Universe))
FALSE  TRUE
13777  1299
> print(table(completetstdf$Group))
Attenuated   Outgroup
       108       1014

## check unigroupdf
> print(table(unigroupdf$universe))
FALSE  TRUE
13095  1981
> print(table(unigroupdf$group))
Attenuated   Outgroup
      1096        885

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
message("Filtering results")
start_time <- Sys.time()
unigroupdf <- universegroup(completedf, expdf, filterdf)
end_time <- Sys.time()
message("\t\t ## Analysis performed in: ", end_time - start_time) # nolint
if (!testonerep) {
  saveRDS(unigroupdf, "/g/romebioinfo/tmp/downstream/unigroupdf.rds")
} else {
  saveRDS(unigroupdf, "/g/romebioinfo/tmp/downstream/unigroupdf-onerep.rds")
}
