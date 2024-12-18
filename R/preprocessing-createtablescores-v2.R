createtablescores <- function(tmpfold, exptabpath, verbose) {

    ## Reading the information about experiments
    if (verbose) message("Reading the information about experiments")
    exptab <- read.csv(exptabpath, header = TRUE)

    ## Retrieving the file paths
    filevec <- list.files(tmpfold, full.names = TRUE)

    ## Splitting the files by experiment names
    expnamevec <- sapply(strsplit(basename(filevec), "-"), "[", 1)

    if (!isTRUE(all.equal(length(unique(table(expnamevec))), 1)))
        stop("Experiments have a different number of files. This should not",
            "happen. Contact the developer.")

    explist <- split(filevec, factor(expnamevec))

    ## Merging files by experiment and direction
    if (verbose) message("\t Merging files by experiment and direction")
    rowidreslist <- mapply(function(currentfiles, currentname, tmpfold,
        verbose) {
            destfile <- file.path(tmpfold, paste0(currentname, ".tsv"))
            if (verbose) message("\t\t Combining all ", currentname,
                " files into ", destfile)
            cmd <- paste0("cat ", paste(currentfiles, collapse = " "), " > ",
                destfile)
            system(cmd)
            return(destfile)
        }, explist, names(explist), MoreArgs = list(tmpfold, verbose))

!!
    colnamevec <- c("biotype", "chr", "coor1", "coor2", "transcript", "gene",
        "strand", "window", "id", "dataset", "score")
    test <- read.delim(rowidreslist[[1]], header = FALSE, sep = "\t",
        na.strings = "NA", dec = ".", col.names = colnamevec,
        stringsAsFactors = FALSE)

    !!start join command - full join is too difficult/long in bash. reuse the purrr::reduce(rowidreslist, dplyr::full_join,
    !! sort properly the final table
}