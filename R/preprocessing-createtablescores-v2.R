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