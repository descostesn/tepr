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

#' Join Bedgraph Files for Protein-Coding and lncRNA Data
#'
#' @description
#' The `joinfiles` function processes bedgraph files located in the specified
#' working directory, joins data related to protein-coding and lncRNA
#' annotations, and outputs a combined result as a TSV file. The function
#' retrieves bedgraph files matching a specific pattern, processes the files to
#' create windows-based summaries, and merges annotations for protein-coding and
#' lncRNA biotypes. The resulting data is written to an output file.
#'
#' @usage
#' joinfiles(workingdir = ".", window = 200, bgpattern = "*.bg",
#'   protscoredir = "protein_coding_score", lncscoredir = "lncRNA_score",
#'   outtsv = "dTAG_Cugusi_stranded_20230810.tsv", verbose = TRUE)
#'
#' @param workingdir The directory containing bedgraph files. Defaults to the
#'  current working directory (`"."`).
#' @param window The window size used for joining the score files. Defaults to
#'  200.
#' @param bgpattern A file pattern to identify bedgraph files. Defaults to
#'  `"*.bg"`.
#' @param protscoredir Directory containing the protein-coding score files.
#'  Defaults to `"protein_coding_score"`.
#' @param lncscoredir Directory containing the lncRNA score files. Defaults to
#'  `"lncRNA_score"`.
#' @param outtsv The output TSV filename where the merged data will be saved.
#'  Defaults to `"dTAG_Cugusi_stranded_20230810.tsv"`.
#' @param nbcpu An integer specifying the number of CPU cores to use for
#'    parallel computation. Default is \code{1}.
#' @param verbose Logical flag to enable verbose output during the function
#'  execution. Defaults to `TRUE`.
#'
#' @return The data.frame with the complete set of annotations and scores.
#'
#' @importFrom purrr map reduce
#' @importFrom tools file_path_sans_ext
#' @importFrom dplyr left_join bind_rows filter
#' @importFrom utils read.delim write.table
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#'   joinfiles(workingdir = "data", window = 100, bgpattern = "*.bedgraph",
#'     protscoredir = "prot_scores", lncscoredir = "lnc_scores",
#'     outtsv = "results.tsv")
#' }
#'
#' @export

joinfiles <- function(workingdir = ".", window = 200, bgpattern = "*.bg", # nolint
    protscoredir = "protein_coding_score", lncscoredir = "lncRNA_score",
    outtsv = "dTAG_Cugusi_stranded_20230810.tsv", nbcpu = 1, verbose = TRUE) {

        ## Defining vectors for protein-coding and lncRNA files
        scoredirvec <- c(protscoredir, lncscoredir)

        ## Retrieving all bedgraph files
        if (verbose) message("Retrieving all bedgraph file paths")
        bedgraphfiles <- list.files(workingdir, pattern = bgpattern,
            full.names = TRUE)

        ## Joining protein coding and lncRNA
        if (verbose) message("Joining protein coding and lncRNA")
        joineddflist <- lapply(scoredirvec, function(scoredir, bedgraphfiles,
            window, nbcpu) {

                if (verbose) message("\t processing ", scoredir)

                files <- bedgraphfiles %>% purrr::map(~{
                    filename <- tools::file_path_sans_ext(basename(.))
                    file.path(scoredir, paste0(filename, ".window", window,
                        ".MANE.wmean.name.score"))})

                ## Reading all files
                if (verbose) message("\t\t Reading all files")
                colnamevec <- c("biotype", "chr", "coor1", "coor2",
                    "transcript", "gene", "strand", "window", "id",
                    "dataset", "score")
                    dflist <- parallel::mclapply(files, read.delim,
                        header = FALSE, sep = "\t", na.strings = "NAN",
                        dec = ".", col.names = colnamevec,
                        stringsAsFactors = FALSE, mc.cores = nbcpu)

                ## Joining all files
                if (verbose) message("\t\t Joining all files")
                joincolvec <- c("biotype", "chr", "coor1", "coor2",
                    "transcript", "gene", "strand", "window", "id")
                ## the last filter remove the PAR genes (pseudoautosomal genes
                ## both in X and Y)
                joineddf <- purrr::reduce(dflist, dplyr::left_join,
                    by = joincolvec) %>% dplyr::filter(strand != "Y")

                return(joineddf)

            }, bedgraphfiles, window, nbcpu)

            ## Merging protein-coding and lncRNA annotations
            if (verbose) message("Merging protein-coding and lncRNA ",
                "annotations")
            bounddf <- purrr::reduce(joineddflist, dplyr::bind_rows)

            ## Writting result table
            outfile <- file.path(workingdir, outtsv)
            if (verbose) message("Writing the result table to ", outfile)
            write.table(bounddf, file = outfile, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
            return(bounddf)
}
