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

createtablescores <- function(bedgraphlistwmean, nbcpubg, exptabpath,
    saveobjectpath = NA, verbose = TRUE) {

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- read.csv(exptabpath, header = TRUE)

        ## Creating a rowid that will be used for merging
        if (verbose) message("\t Adding rowid for each bedgraph")
        rowidreslist <- .createrowidlist(bedgraphlistwmean, nbcpubg)

        if (verbose) message("\t Joining the elements of each bedgraph")
        df <- purrr::reduce(rowidreslist, dplyr::full_join,
            by = c("chrom", "start.window", "end.window", "strand.window",
                "gene", "biotype.window", "window", "coord", "transcript",
                "rowid"))

!!!!!!! 
        if (verbose) message("\t Sorting columns")
        test <- df %>% relocate(c("biotype.window", "chrom", "start.window",
        "end.window", "transcript", "gene", "strand.window", "window", "rowid"))
        test <- test[, -which(colnames(test) == "coord")]


- NAMES OF EXP

    idxcolscores <- grep("_score", colnames(completeframedf))
    newscorenames <- unlist(apply(exptab, 1, function(x) {
        return(paste0(x["condition"], "_rep", x["replicate"], ".", x["strand"]))
    }, simplify = FALSE))
    colnames(completeframedf)[idxcolscores] <- paste(newscorenames, "score", sep = "_")

- NEW DF OF EXP NAMES
    dfexpnameslist <- lapply(newscorenames, rep, nrow(completeframedf))
    dfexpnames <- do.call("cbind", dfexpnameslist)

- 

- ASSOCIATING SCORES WITH THE COL OF EXP NAME
expnamescorelist <- parallel::mcmapply(function(expnamecol, scorecol) {
    return(cbind(expnamecol, scorecol))
}, dfexpnames, completeframedf[idxcolscores], SIMPLIFY = FALSE,
    mc.silent = FALSE, mc.cores = nbcpubg)
expnamescoredf <- do.call("cbind", expnamescorelist)

-  COMBINING THE TWO FINAL DF
finaldf <- cbind(noexpdf, expnamescoredf) 
!!!!!!!!!!!!!

        if (!is.na(saveobjectpath)) {
            outfile <- file.path(saveobjectpath, "finaltab.rds")
            if (verbose) message("\t Saving ", outfile)
            saveRDS(completeframedf, file = outfile)
        }

        return(completeframedf)
}
