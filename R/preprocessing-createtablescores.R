# .createrowidlist <- function(bedgraphlistwmean, nbcpubg) { # nolint

#         rowidreslist <- parallel::mclapply(bedgraphlistwmean, function(tab) {
#             ## Create rowid string
#             rowidvec <- paste(tab$transcript, tab$gene, tab$strand, tab$window,
#                 sep = "_")
#             ## Inserting rowid col after window
#             tab <- tab %>% tibble::add_column(rowid = rowidvec,
#                 .after = "window")
#             return(tab)
#         }, mc.cores = nbcpubg)

#     return(rowidreslist)
# }

.orderingtable <- function(df, exptab, verbose) { # nolint

    # if (verbose) message("\t\t Sorting and renaming information columns")
    # df <- df %>% dplyr::relocate(biotype, .before = chrom) # nolint
    # idxtorename <- match(c("chrom", "start", "end", "rowid"), colnames(df))
    # colnames(df)[idxtorename] <- c("chr", "coor1", "coor2", "id")

    if (verbose) message("\t\t Renaming score columns")
    # idxcolscores <- grep("_score", colnames(df))
    # expcolnames <- unlist(apply(exptab, 1, function(x) {
    #     return(paste0(x["condition"], "_rep", x["replicate"], ".",
    #         x["strand"]))
    # }, simplify = FALSE))
    # newscorenames <- paste(expcolnames, "score", sep = "_")
    # colnames(df)[idxcolscores] <- newscorenames

    # if (verbose) message("\t\t Creating experiment columns")
    # ## The format of the experiment column is title "HS_rep1.plus", content "HS_rep1.forward" # nolint
    # directionexpstr <- unlist(apply(exptab, 1, function(x) {
    #     return(paste0(x["condition"], "_rep", x["replicate"], ".",
    #         x["direction"]))}, simplify = FALSE))
    # dfexpnameslist <- lapply(directionexpstr, rep, nrow(df))
    # dfexpnames <- do.call("cbind", dfexpnameslist)
    # colnames(dfexpnames) <- expcolnames

    # if (verbose) message("\t\t Combining the experiment cols to the table")
    # df <- cbind(df, dfexpnames)
    # df <- tibble::as_tibble(df)

    # if (verbose) message("\t\t Placing exp name columns before corresponding",
    #     " scores")
    # df <- purrr::reduce(.x = purrr::map2(expcolnames, newscorenames, c),
    #     .f = ~ dplyr::relocate(.x, .y[1], .before = .y[2]), .init = df)

    return(df)
}


#' Create a Table of Scores from Bedgraph Data
#'
#' @description
#' This function processes the scores of bedgraph files retrieved on transcripts
#' to create a merged table of scores. It performs tasks such as adding unique
#' row IDs, joining bedgraph data, reordering and renaming columns, and
#' preparing a final table for further analysis.
#'
#' @usage
#' createtablescores(bedgraphlistwmean, nbcpubg, exptabpath,
#'  saveobjectpath = NA, reload = FALSE, showtime = TRUE, verbose = TRUE)
#'
#' @param bedgraphlistwmean A list for each bedgraph with values on transcript
#'  windows. It was obtained with the function blacklisthighmap.
#' @param nbcpubg Number of CPU cores to use for bedgraph-level operations.
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param saveobjectpath Path to save intermediate R objects. Default is `NA`
#'  and R objects are not saved.
#' @param reload Logical. If `TRUE`, reloads existing saved objects to avoid
#'  recomputation. Default is `FALSE`. If the function failed during object
#'  saving, make sure to delete the corresponding object.
#' @param showtime Logical. If `TRUE`, displays timing information. Default is
#'  `FALSE`.
#' @param verbose Logical. If `TRUE`, provides detailed messages during
#'  execution. Default is `TRUE`.
#'
#' @return A tibble containing the final merged and processed table with scores
#' and experiment details.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads experiment details from the provided exptabpath.
#' 2. Adds a unique row ID to each bedgraph data frame based on transcript,
#'  gene, strand, and window.
#' 3. Merges all bedgraph data frames into a single data frame using the unique
#'  row ID.
#' 4. Reorders and renames columns for consistency and clarity.
#' 5. Adds experiment-related columns based on the provided experiment table.
#' 6. Optionally saves the resulting table to the specified path.
#'
#' @examples
#' # Example usage of createtablescores function:
#' bedgraphlist <- list(bedgraph1, bedgraph2, bedgraph3)
#' exptabpath <- "path/to/experiment_table.csv"
#' output_path <- "path/to/save_directory"
#' final_table <- createtablescores(bedgraphlist, nbcpubg = 4, exptabpath,
#'  saveobjectpath = output_path, reload = FALSE, showtime = TRUE,
#'  verbose = TRUE)
#'
#' @seealso
#' [blacklisthighmap]
#'
#' @importFrom parallel mclapply
#' @importFrom purrr reduce map2
#' @importFrom dplyr full_join relocate
#' @importFrom tibble add_column as_tibble
#'
#' @export

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

createtablescores <- function(tmpfold, verbose) {

    ## Retrieving the file paths
    filevec <- list.files(tmpfold, full.names = TRUE)

    ## Splitting the files by experiment names
    expnamevec <- sapply(strsplit(basename(filevec), "-"), "[", 1)

    if (!isTRUE(all.equal(length(unique(table(expnamevec))), 1)))
        stop("Experiments have a different number of files. This should not",
            "happen. Contact the developer.")

    explist <- split(filevec, factor(expnamevec))

    if(verbose) message("\t Merging files by experiment and direction")
    !! check if the direction is differentiated
    rowidreslist <- mapply(function(currentfiles, currentname, tmpfold,
        verbose) {
            destfile <- file.path(tmpfold, paste0(currentname, ".tsv"))
            if (verbose) message("\t\t Combining all ", currentname,
                " files into ", destfile)
            cmd <- paste0("cat ", paste(currentfiles, collapse = " "), " > ",
                destfile)
        #   !! TO UNCOMMENT system(cmd)
        !! reading the merged table and returning it - see how much memory it takes

            return(destfile)
        }, explist, names(explist), MoreArgs = list(tmpfold, verbose))
    
     
    !!start join command - full join is too difficult/long in bash. reuse the purrr::reduce(rowidreslist, dplyr::full_join,
    !! sort properly the final table
}