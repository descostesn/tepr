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
