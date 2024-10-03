.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname, # nolint
                " should be unique, contact the developer.") # nolint
}

.extractstr <- function(transtable) {

    str <- as.character(unique(transtable$strand))
    .checkunique(str, "str")
    if (isTRUE(all.equal(str, "+"))) {
        str <- "plus"
    } else if (isTRUE(all.equal(str, "-"))) {
        str <- "minus"
    } else {
        stop("In .computeecdf, strand is neither plus or minus. This ",
            "should not happen. Contact the developer.")
    }
    return(str)
}

.colnamecheck <- function(colnamevec, tab) {
            invisible(sapply(colnamevec, function(currentcol, tab) {
            idx <- grep(currentcol, colnames(tab))
            if (isTRUE(all.equal(length(idx), 0)))
                stop("The column ", currentcol, " does not exist in the ",
                    "provided table.")
        }, tab))
}





joinfiles <- function(workingdir = ".", window = 200, bgpattern = "*.bg",
    protscoredir = "protein_coding_score", lncscoredir = "lncRNA_score",
    outtsv = "dTAG_Cugusi_stranded_20230810.tsv", verbose = TRUE) {

        ## Defining vectors for protein-coding and lncRNA files
        scoredirvec <- c(protscoredir, lncscoredir)

        ## Retrieving all bedgraph files
        if (verbose) message("Retrieving all bedgraph file paths")
        bedgraphfiles <- list.files(workingdir, pattern = bgpattern,
            full.names = TRUE)

        ## Joining protein coding and lncRNA
        if (verbose) message("Joining protein coding and lncRNA")
        joineddflist <- lapply(scoredirvec, function(scoredir, bedgraphfiles,
            window) {

                if (verbose) message("\t processing ", scoredir)

                files <- bedgraphfiles %>% purrr::map(~{
                    filename <- tools::file_path_sans_ext(basename(.))
                    file.path(scoredir, paste0(filename, ".window", window,
                        ".MANE.wmean.name.score"))})

                # Reading all the files
                if (verbose)
    colnamevec <- c("biotype", "chr", "coor1", "coor2", "transcript", "gene",
        "strand", "window", "id", "dataset", "score")
    dflist <- lapply(files, read.delim, header = FALSE, sep = "\t",
        na.strings = "NAN", dec = ".", col.names = colnamevec,
        stringsAsFactors = FALSE)

    # joining all the files
    joincolvec <- c("biotype", "chr", "coor1", "coor2", "transcript", "gene",
        "strand", "window", "id")
    ## the last filter remove the PAR genes (pseudoautosomal genes both in X and Y)
    joineddf <- purrr::reduce(dflist, dplyr::left_join, by = joincolvec) %>%
        dplyr::filter(strand != "Y")

    return(joineddf)

}, bedgraphfiles, window)


}


# joining all the files
bounddf <- purrr::reduce(joineddflist, dplyr::bind_rows)

write.table(bounddf, file = file.path(workingdir, outtsv), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
