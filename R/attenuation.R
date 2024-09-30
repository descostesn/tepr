.summarytrans <- function(bytransmeanlist, nbcpu) {
  summarydflist <- mclapply(bytranslistmean, function(trans) {
    coor1 <- min(trans$coor1)
    coor2 <- max(trans$coor2)
    return(data.frame(chr = trans$chr[1], coor1, coor2,
          strand = trans$strand[1], gene = trans$gene[1],
          transcript = trans$transcript[1], size = coor2 - coor1 + 1))
  }, mc.cores = nbcpu)
  summarydf <- do.call("rbind", summarydflist)
  return(summarydf)
}

.computeupdown <- function(completbytrans, condvec, nbcpu) {

  updownbytranslist <- mclapply(completbytrans, function(trans, condvec) {

    ## Ordering by coordinates (security)
    trans <- trans[order(trans$coord), ]

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
      colnames(res) <- c("transcript", paste0("UP_mean_", cond),
              paste0("DOWN_mean_", cond), paste0("Attenuation_", cond))
      return(res)
    }, trans)

    return(do.call("cbind", updownlist))
  }, condvec, mc.cores = nbcpu)

  updowndf <- do.call("rbind", updownbytranslist)
  updowndf <- updowndf[, -which(duplicated(colnames(updowndf)))]
  return(updowndf)
}

.filterattenuation <- function(auckneenasumatt, condvec, pval, replaceval, # nolint
    verbose) {

        mat <- auckneenasumatt
        if (verbose) message("\t\t Replacing non-significant attenuations by ",
            replaceval)
        invisible(sapply(condvec, function(cond, replaceval) {
            pauccond <- paste0("p_AUC_", cond)
            ## Replacing Attenuation value if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("Attenuation_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("Attenuation_", cond)]]))
            ## Replacing knee values if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("knee_AUC_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("knee_AUC_", cond)]]))
        }, replaceval))

        return(mat)
}

attenuation <- function(allaucdf, kneedf, matnatrans, bytranslistmean, expdf,
  dfmeandiff, nbcpu = 1, significant = FALSE, replaceval = NA, pval = 0.1,
  verbose = TRUE) {

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
      if (verbose) message("\t Merging attenuation to the complete table")
      auckneenasumatt <- merge(auckneenasum, updowndf, by = "transcript")

      ## Replace the attenuation values by replaceval if p_AUC_cond >= pval
      if (significant) {
        if (verbose) message("\t Keeping significant attenuation")
        auckneenasumatt <- .filterattenuation(auckneenasumatt, condvec, pval,
            replaceval, verbose)
      }

      return(auckneenasumatt)
}
