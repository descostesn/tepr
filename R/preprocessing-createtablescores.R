.mergefiles <- function(explist, tmpfold, verbose) {

    if (verbose) message("Merging files by experiment and direction")
    mergedfilelist <- mapply(function(currentfiles, currentname, tmpfold,
        verbose) {
            destfile <- file.path(tmpfold, paste0(currentname, ".tsv"))
            if (verbose) message("\t Combining all ", currentname,
                " files into ", destfile)
            cmd <- paste0("cat ", paste(currentfiles, collapse = " "),
                " > ", destfile)
            system(cmd)
            return(destfile)
        }, explist, names(explist), MoreArgs = list(tmpfold, verbose))

    return(mergedfilelist)
}

.fulljoinfun <- function(mergedfilelist, idxvec, showmemory, verbose) {

    colnamevec <- c("biotype", "chr", "coor1", "coor2", "transcript", "gene",
        "strand", "window", "id", "dataset", "score")
    colnamejoin <- colnamevec[-c(10, 11)] ## Remove dataset and score
    if (verbose) message("Reading files and joining to the final table")
    firstpath <- mergedfilelist[[idxvec[1]]]
    if (verbose) message("\t Reading and joining ", firstpath)
    finaltab <- read.delim(firstpath, header = FALSE,
        sep = "\t", na.strings = "NA", dec = ".", col.names = colnamevec,
        stringsAsFactors = FALSE)

    for (idx in idxvec[-1]) {
        currentpath <- mergedfilelist[[idx]]
        if (verbose) message("\t Reading and joining ", currentpath)
        tab <- read.delim(currentpath, header = FALSE, sep = "\t",
            na.strings = "NA", dec = ".", col.names = colnamevec,
            stringsAsFactors = FALSE)
        finaltab <- dplyr::full_join(finaltab, tab, by = colnamejoin)
        rm(tab)
        if (showmemory) print(gc()) else invisible(gc())
    }
    return(finaltab)
}

#' Create a Unified Table of Scores
#'
#' @description
#' This function processes and combines table scores of each bedgraph and each
#' chromosome stored in the temporary folder into a unified table.
#'
#' @usage
#' createtablescores(tmpfold, exptabpath, showmemory = FALSE, showtime = TRUE,
#'   savefinaltable = TRUE, finaltabpath = "./", finaltabname = "anno.tsv",
#'  verbose)
#'
#' @param tmpfold A string specifying the temporary folder containing the
#'  score files created with the function 'blacklisthighmap'.
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param showmemory Logical; if `TRUE`, memory usage is printed during
#'  processing. Default is `FALSE`.
#' @param showtime Logical; if `TRUE`, the execution time of the function is
#'  printed. Default is `TRUE`.
#' @param savefinaltable Logical; if `TRUE`, the resulting table is saved to
#'  disk. Default is `TRUE`.
#' @param finaltabpath A string specifying the directory where the final table
#'  should be saved. Default is `"./"`.
#' @param finaltabname A string specifying the name of the final table file.
#'  Default is `"anno.tsv"`.
#' @param verbose Logical; if `TRUE`, detailed messages are printed during
#'  execution.
#'
#' @return A data frame containing the unified table of scores.
#'
#' @details
#' This function first merges files belonging to the same experiment and
#' direction. These files are combined into a single table providing two columns
#' per experiment. The first gives the name of the experiment and the second the
#' scores. The resulting table also includes annotations for each transcript.
#'
#' @examples
#' # Example usage:
#' tmpfold <- "path/to/temp/folder"
#' exptabpath <- "path/to/experiment_table.csv"
#' finaltab <- createtablescores(tmpfold = tmpfold, exptabpath = exptabpath,
#'   showmemory = TRUE, showtime = TRUE, savefinaltable = TRUE,
#'   finaltabpath = "./results", finaltabname = "final_scores.tsv",
#'   verbose = TRUE)
#'
#' @importFrom dplyr full_join
#'
#' @seealso
#' [blacklisthighmap]
#'
#' @export

createtablescores <- function(tmpfold, exptabpath, showmemory = FALSE,
    showtime = TRUE, savefinaltable = TRUE, finaltabpath = "./",
    finaltabname = "anno.tsv", verbose = TRUE) {

        if (showtime) start_time_fun <- Sys.time()
        if (verbose) message("\n ## Merging results of each bedgraph into a ",
            "single table ##\n")

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- read.csv(exptabpath, header = TRUE)

        ## Retrieving the file paths
        filevec <- list.files(tmpfold, full.names = TRUE)

        ## Splitting the files by experiment names
        expnamevec <- sapply(strsplit(basename(filevec), "-"), "[", 1)

        if (!isTRUE(all.equal(length(unique(table(expnamevec))), 1)))
            stop("Experiments have a different number of files. This should ",
                "not happen. Contact the developer.")

        explist <- split(filevec, factor(expnamevec))

        ## Merging files by experiment and direction
        mergedfilelist <- .mergefiles(explist, tmpfold, verbose)

        ## Retrieving the exp name in the right order from exptab
        orderedexpvec <- paste0(exptab$condition, exptab$replicate,
            exptab$direction)
        idxvec <- match(orderedexpvec, names(mergedfilelist))
        idxna <- which(is.na(idxvec))
        if (!isTRUE(all.equal(length(idxna), 0)))
            stop("The merged file names do not correspond to the exptab. This",
                "should not happen. Contact the developer.")

        ## Reading each merged file and combining it to the final table
        finaltab <- .fulljoinfun(mergedfilelist, idxvec, showmemory, verbose)

        ## Filling the experiment columns
        if (verbose) message("Filling the experiment columns")
        idxdatasetvec <- grep("dataset", colnames(finaltab))
        nbrows <- nrow(finaltab)
        for (idxdataset in idxdatasetvec) {
            expname <- as.character(na.omit(unique(finaltab[, idxdataset])))
            finaltab[, idxdataset] <- rep(expname, nbrows)
        }

        if (savefinaltable) {
            if (!file.exists(finaltabpath))
                dir.create(finaltabpath, recursive = TRUE)
            outfile <- file.path(finaltabpath, finaltabname)
            if (verbose) message("\n ## Saving the final table to ", outfile)
            write.table(finaltab, file = outfile, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
        }

        if (showtime) {
            end_time_fun <- Sys.time()
            timing <- end_time_fun - start_time_fun
            message("\t\t ## Final table created in: ",
                format(timing, digits = 2))
        }

        return(finaltab)
}
