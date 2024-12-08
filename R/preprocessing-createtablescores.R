.createrowidlist <- function(bedgraphlistwmean, nbcpubg) { # nolint

        rowidreslist <- parallel::mclapply(bedgraphlistwmean, function(tab) {

        rowidvec <- paste(tab$transcript, tab$gene, tab$strand, tab$window,
            sep = "_")
        # ENST00000000233.10_ARF5_+_1
        # protein-coding_chr7_127588411_127588427_+_ARF5_ENST00000000233.10_frame1_coord1 # nolint
        # TO REMOVE rowidvec <- paste(tab$biotype.window, tab$chrom, tab$start.window, tab$end.window, tab$strand.window, tab$gene, tab$transcript, paste0("frame", tab$window), paste0("coord", tab$coord), sep = "_") # nolint

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
!!!!!!!!!!!!!!!!!!! TO CONTINUE
        tab <- dplyr::relocate(tab, colnamevec[idxscore], .after = "coord")
        return(tab)
    }, mc.cores = nbcpubg)

    return(rowidreslist)
}

.orderingtable <- function(df, exptab, verbose) {

    if (verbose) message("\t\t Sorting columns")
    orderedcolvec <- c("biotype.window", "chrom", "start.window",
        "end.window", "transcript", "gene", "strand.window", "window",
        "rowid")
    df <- df %>% dplyr::relocate(orderedcolvec)
    df <- df[, -which(colnames(df) == "coord")]

    if (verbose) message("\t\t Renaming score columns")
    idxcolscores <- grep("_score", colnames(df))
    expcolnames <- unlist(apply(exptab, 1, function(x) {
        return(paste0(x["condition"], "_rep", x["replicate"], ".",
            x["strand"]))
    }, simplify = FALSE))
    newscorenames <- paste(expcolnames, "score", sep = "_")
    colnames(df)[idxcolscores] <- newscorenames

    if (verbose) message("\t\t Creating experiment columns")
    dfexpnameslist <- lapply(expcolnames, rep, nrow(df))
    dfexpnames <- do.call("cbind", dfexpnameslist)
    colnames(dfexpnames) <- expcolnames

    if (verbose) message("\t\t Combining the experiment cols to the table")
    df <- cbind(df, dfexpnames)
    df <- tibble::as_tibble(df)

    if (verbose) message("\t\t Placing exp name columns before corresponding",
        " scores")
    df <- purrr::reduce(.x = purrr::map2(expcolnames, newscorenames, c),
        .f = ~ dplyr::relocate(.x, .y[1], .before = .y[2]), .init = df)

    return(df)
}

createtablescores <- function(bedgraphlistwmean, nbcpubg, exptabpath,
    saveobjectpath = NA, reload = FALSE, showtime = TRUE, verbose = TRUE) {

        if (showtime) start_time_fun <- Sys.time()
        dfobj <- file.path(saveobjectpath, "finaltab.rds")

        if (!reload || !file.exists(dfobj)) {

            ## Reading the information about experiments
            if (verbose) message("Reading the information about experiments")
            exptab <- read.csv(exptabpath, header = TRUE)

            ## Creating a rowid that will be used for merging
            if (verbose) message("\t Adding rowid for each bedgraph")
            rowidreslist <- .createrowidlist(bedgraphlistwmean, nbcpubg)

            if (verbose) message("\t Joining the elements of each bedgraph")
            df <- purrr::reduce(rowidreslist, dplyr::full_join,
                by = c("chrom", "start.window", "end.window", "strand.window",
                    "gene", "biotype.window", "window", "transcript", "rowid"))
                # TO REMOVE by = c("chrom", "start.window", "end.window", "strand.window", "gene", "biotype.window", "window", "coord", "transcript", "rowid")) # nolint

            rm(rowidreslist)
            invisible(gc())

            if (verbose) message("\t Preparing final table")
            df <- .orderingtable(df, exptab, verbose)

            if (!is.na(saveobjectpath)) {
                if (verbose) message("\t Saving ", dfobj)
                saveRDS(df, file = dfobj)
            }
        } else {

            if (verbose) message("Loading object ", dfobj)
            df <- readRDS(dfobj)
        }

        if (showtime) {
            end_time_fun <- Sys.time()
            timing <- end_time_fun - start_time_fun
            message("\t\t ## Final table created in: ", timing) # nolint
        }

        return(df)
}
