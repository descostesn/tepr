.createrowidlist <- function(bedgraphlistwmean, nbcpubg) { # nolint

        rowidreslist <- parallel::mclapply(bedgraphlistwmean, function(tab) {
            ## Create rowid string
            rowidvec <- paste(tab$transcript, tab$gene, tab$strand, tab$window,
                sep = "_")
            ## Inserting rowid col after window
            tab <- tab %>% tibble::add_column(rowid = rowidvec,
                .after = "window")
            return(tab)
        }, mc.cores = nbcpubg)

    return(rowidreslist)
}

.orderingtable <- function(df, exptab, verbose) {

!!!!!!!!!!!!!!!!!
    nic:
    [1] "chrom"              "biotype"            "start"
 [4] "end"                "transcript"         "gene"
 [7] "strand"             "window"             "rowid"
[10] "ctrl1forward_score" "ctrl1reverse_score" "ctrl2forward_score"
[13] "ctrl2reverse_score" "HS1forward_score"   "HS1reverse_score"
[16] "HS2forward_score"   "HS2reverse_score"
    vic:
    [1] "biotype"               "chr"                   "coor1"
 [4] "coor2"                 "transcript"            "gene"
 [7] "strand"                "window"                "id"
[10] "ctrl_rep1.plus"        "ctrl_rep1.plus_score"  "ctrl_rep1.minus"
[13] "ctrl_rep1.minus_score" "ctrl_rep2.plus"        "ctrl_rep2.plus_score"
[16] "ctrl_rep2.minus"       "ctrl_rep2.minus_score" "HS_rep1.plus"
[19] "HS_rep1.plus_score"    "HS_rep1.minus"         "HS_rep1.minus_score"
[22] "HS_rep2.plus"          "HS_rep2.plus_score"    "HS_rep2.minus"
[25] "HS_rep2.minus_score"
    !!!!!!!!!!!!!!!!

    if (verbose) message("\t\t Sorting and renaming information columns")
    df <- df %>% dplyr::relocate(biotype, .before = chrom)
    idxtorename <- match(c("chrom", "start", "end", "rowid"), colnames(df))
    colnames(df)[idxtorename] <- c("chr", "coor1", "coor2", "id")

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

            ## Joining the elements of each bedgraph
            if (verbose) message("\t Joining the elements of each bedgraph")
            if (showtime) start_time_join <- Sys.time()
            df <- purrr::reduce(rowidreslist, dplyr::full_join,
                by = c("biotype", "chrom", "start", "end", "transcript", "gene",
                    "strand", "window", "rowid"))
            if (showtime) {
                end_time_join <- Sys.time()
                timing <- end_time_join - start_time_join
            message("\t\t ## Joined table in: ", timing) # nolint
            }

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
