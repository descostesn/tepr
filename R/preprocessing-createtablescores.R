.createrowidlist <- function(bedgraphlistwmean, nbcpubg) { # nolint

        rowidreslist <- parallel::mclapply(bedgraphlistwmean, function(tab) {

        rowidvec <- paste(tab$biotype.window, tab$chrom, tab$start.window,
            tab$end.window, tab$strand.window, tab$gene,
            tab$transcript, paste0("frame", tab$window),
            paste0("coord", tab$coord), sep = "_")

        ## Inserting rowid col after transcript
        tab <- tab %>% tibble::add_column(rowid = rowidvec,
            .after = "transcript")

        ## Remove bedgraph columns
        tab <- tab[, -grep(".bg", colnames(tab))]

        ## Move the score column at the end of the table
        colnamevec <- colnames(tab)
        idxscore <- grep("_score", colnamevec)
        if (!isTRUE(all.equal(length(idxscore), 1)))
            stop("When creating the final table, the score column is not ",
                "unique for a given bedgraph. This should not happen. Contact",
                " the developer.")
        tab <- dplyr::relocate(tab, colnamevec[idxscore], .after = "coord")
        return(tab)
    }, mc.cores = nbcpubg)

    return(rowidreslist)
}

createtablescores <- function(bedgraphlistwmean, nbcpubg, saveobjectpath = NA, # nolint
    verbose = TRUE) {

    ## Creating a rowid that will be used for merging
    if (verbose) message("\t Adding rowid for each bedgraph")
    rowidreslist <- .createrowidlist(bedgraphlistwmean, nbcpubg)

    if (verbose) message("\t Joining the elements of each bedgraph")
    completeframedf <- purrr::reduce(rowidreslist, dplyr::full_join,
        by = c("chrom", "start.window", "end.window", "strand.window", "gene",
        "biotype.window", "window", "coord", "transcript", "rowid"))

    if (!is.na(saveobjectpath))
        saveRDS(completeframedf, file = file.path(saveobjectpath,
            "finaltab.rds"))

    return(completeframedf)
}
