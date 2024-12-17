.createallannobed <- function(exptabpath, gencodepath, saveobjectpath, reload,
    showtime, verbose) {

        if (verbose) message("## Filtering gencode annotations ##\n")
        allannobedobjpath <- file.path(saveobjectpath, "allannobed.rds")
        if (!reload || !file.exists(allannobedobjpath)) {
            allannobed <- retrieveanno(exptabpath, gencodepath, saveobjectpath,
                showtime, verbose)
        } else {
            if (verbose) message("Loading ", allannobedobjpath)
            allannobed <- readRDS(allannobedobjpath)
        }
        return(allannobed)
}

.createallwindowsbed <- function(allannobed, windsize, nbcputrans, showtime,
    saveobjectpath, reload, verbose) {

        if (verbose) message("\n ## Splitting transcripts into windows ##\n")
            allwindowsbedobjpath <- file.path(saveobjectpath,
                "allwindowsbed.rds")
        if (!reload || !file.exists(allwindowsbedobjpath)) {
            allwindowsbed <- makewindows(allannobed, windsize, nbcputrans,
                verbose, saveobjectpath, showtime)
        } else {
            if (verbose) message("Loading ", allwindowsbedobjpath)
            allwindowsbed <- readRDS(allwindowsbedobjpath)
        }

        return(allwindowsbed)
}

.createbedgraphlistwmean <- function(maptrackpath, blacklistshpath,
    exptabpath, nbcputrans, allwindowsbed, windsize, showtime, saveobjectpath,
    reload, verbose) {

        if (verbose) message("\n ## Retrieving the values of the bedgraph ",
            "files, removing black lists and keeping scores landing on high ",
            "mappability intervals ##\n")
        bedgraphlistwmeanobjpath <- file.path(saveobjectpath,
            "bedgraphlistwmean.rds")

        if (!reload || !file.exists(bedgraphlistwmeanobjpath)) {
            bedgraphlistwmean <- blacklisthighmap(maptrackpath, blacklistshpath,
                exptabpath, nbcputrans, allwindowsbed, windsize, saveobjectpath,
                reload, showtime, verbose)
        } else {
            if (verbose) message("Loading ", bedgraphlistwmeanobjpath)
                bedgraphlistwmean <- readRDS(bedgraphlistwmeanobjpath)
        }
        return(bedgraphlistwmean)
}

#' Preprocess Genomic Data
#'
#' @description
#' This function preprocesses genomic data by filtering, splitting, and
#' annotating transcripts based on the given parameters. It generates a final
#' annotated table.
#'
#' @usage
#' preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
#' blacklistshpath, nbcputrans = 1, nbcpubg = 1, finaltabpath = "./",
#'     finaltabname = "anno.tsv", saveobjectpath = NA, savefinaltable = TRUE,
#'     reload = FALSE, showtime = FALSE, verbose = TRUE)
#'
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param gencodepath Path to the GENCODE annotation file.
#' @param windsize Window size for splitting transcripts into intervals.
#' @param maptrackpath Path to the mappability track file.
#' @param blacklistshpath Path to the blacklist file.
#' @param nbcputrans Number of CPU cores to use for transcript-level operations.
#' @param nbcpubg Number of CPU cores to use for bedgraph-level operations.
#' @param finaltabpath Directory path to save the final annotated table.
#' @param finaltabname Name of the final annotated table file.
#' @param saveobjectpath Path to save intermediate R objects. Default is `NA`
#'  and R objects are not saved.
#' @param savefinaltable Logical. If `TRUE`, saves the final table as TSV to
#'  disk. Default is `TRUE`.
#' @param reload Logical. If `TRUE`, reloads existing saved objects to avoid
#'  recomputation. Default is `FALSE`. If the function failed during object
#'  saving, make sure to delete the corresponding object.
#' @param showtime Logical. If `TRUE`, displays timing information. Default is
#'  `FALSE`.
#' @param verbose Logical. If `TRUE`, provides detailed messages during
#'  execution. Default is `TRUE`.
#'
#' @return A data frame containing the final annotated table.
#'
#' @details
#' The `preprocessing` function performs the following steps:
#'
#' 1. Filters GENCODE annotations to retrieve "transcript" entries and
#'  categorizes them as protein-coding or long non-coding RNA.
#' 2. Splits transcripts into windows of size `windsize`, filtering out
#'  intervals smaller than this size and removing entries with "PAR_Y" in their
#'  Ensembl names.
#' 3. Retrieves values from bedgraph files, applies blacklist filtering, and
#'  retains scores in high-mappability regions.
#' 4. Merges results into a final table and optionally saves it to disk.
#'
#' If `reload` is set to `TRUE` and saved objects exist, the function avoids
#' recomputation by reusing those objects.
#'
#' @examples
#' # Example usage:
#' preprocessing(
#'   exptabpath = "experiment_table.tsv",
#'   gencodepath = "gencode.v37.annotation.gtf",
#'   windsize = 200,
#'   maptrackpath = "mappability.bedgraph",
#'   blacklistshpath = "blacklist.bed",
#'   finaltabpath = "./results",
#'   finaltabname = "annotated_results.tsv",
#'   saveobjectpath = "./intermediate_objects",
#'   savefinaltable = TRUE,
#'   reload = FALSE,
#'   showtime = TRUE,
#'   verbose = TRUE
#' )
#'
#' @importFrom tibble as_tibble add_column
#' @importFrom rtracklayer import.bedGraph
#' @importFrom dplyr relocate filter full_join
#' @importFrom valr bed_intersect
#' @importFrom parallel mclapply makeCluster parLapply stopCluster
#' @importFrom purrr reduce map2
#'
#' @export

preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans = 1, nbcpubg = 1,
    finaltabpath = "./", finaltabname = "anno.tsv", saveobjectpath = NA,
    savefinaltable = TRUE, reload = FALSE, showtime = FALSE, verbose = TRUE) {

    if (reload && file.exists(file.path(saveobjectpath, "finaltable.rds")))
        stop("The final table already exists, set reload = FALSE to create",
            "it again.")

    if (showtime) start_time_preprocessing <- Sys.time()

    ## This function filters gencode annotations to retrieve "transcript". It
    ## then distinguishes transcripts coming from protein coding genes
    ## (MANE_Select) and those coming from long non-coding genes (lncRNA,
    ## Ensembl_canonical). It returns the combination of the two types of
    ## transcripts that are distinguished by the column 'biotype'.
    allannobed <- .createallannobed(exptabpath, gencodepath, saveobjectpath,
        reload, showtime, verbose)

    ## This functions uses the annotations filtered from gencode (allannobed).
    ## It removes any ensembl names containing "PAR_Y", filters out intervals
    ## smaller than windsize and splits each transcript into "windsize" windows.
    allwindowsbed <- .createallwindowsbed(allannobed, windsize, nbcputrans,
        showtime, saveobjectpath, reload, verbose)

    if (verbose) message("\t\t Deleting objects and free memory")
    rm(allannobed)
    invisible(gc())

    ## Retrieving the values of the bedgraph files, removing black lists and
    ## keeping scores landing on high mappability intervals
    bedgraphlistwmean <- .createbedgraphlistwmean(maptrackpath, blacklistshpath,
        exptabpath, nbcputrans, allwindowsbed, windsize, genomename, showtime,
        saveobjectpath, reload, verbose)

    ## Creating the final table from the information retrieved from
    ## blacklisthighmap
    if (verbose) message("\n ## Merging results of each bedgraph into a ",
        "single table ##\n")
    finaltable <- createtablescores(bedgraphlistwmean, nbcpubg, exptabpath,
        saveobjectpath, verbose)

    if (savefinaltable) {
        outfile <- file.path(finaltabpath, finaltabname)
        if (verbose) message("\n ## Saving the final table to ", outfile)
        write.table(finaltable, file = outfile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
    }

    !!!!!!! remove all saved obj if set to true (must be the default)

    if (showtime) {
        end_time_preprocessing <- Sys.time()
        timing <- end_time_preprocessing - start_time_preprocessing
            message("\n\n\t ## Total preprocessing in: ", timing) # nolint
    }

    return(finaltable)
}
