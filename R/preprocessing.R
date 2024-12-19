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

.createbedgraphlistwmean <- function(maptrackpath, blacklistshpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, genomename, showtime, showmemory,
    saveobjectpath, reload, tmpfold, verbose) {

        if (verbose) message("\n ## Retrieving the values of the bedgraph ",
            "files, removing black lists and keeping scores landing on high ",
            "mappability intervals ##\n")
        blacklisthighmap(maptrackpath, blacklistshpath, exptabpath,
            nbcputrans, allwindowsbed, windsize, genomename, saveobjectpath,
            tmpfold, reload, showtime, showmemory, verbose)
}



#' Preprocess Experimental Data for Genomic Analysis
#'
#' This function orchestrates a pipeline for preprocessing genomic data,
#' including filtering annotations, splitting transcripts into windows,
#' retrieving bedgraph values, and generating a final annotated table.
#'
#' @usage
#' preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
#' blacklistshpath, genomename, nbcputrans = 1, finaltabpath = "./",
#' finaltabname = "anno.tsv", tmpfold = "./tmp", saveobjectpath = NA,
#' savefinaltable = TRUE, reload = FALSE, showtime = FALSE, showmemory = FALSE,
#' deletetmp = TRUE, verbose = TRUE)
#'
#' @param exptabpath Character. Path to the experiment table file.
#' @param gencodepath Character. Path to the Gencode annotation file.
#' @param windsize Integer. Window size for splitting transcripts.
#' @param maptrackpath Character. Path to the mappability track file.
#' @param blacklistshpath Character. Path to the blacklist file.
#' @param genomename Character. Name of the genome assembly (e.g., "hg38").
#' @param nbcputrans Integer. Number of CPUs to use for transcript processing.
#'  Default is 1.
#' @param finaltabpath Character. Path where the final annotated table will be
#'  saved. Default is "./".
#' @param finaltabname Character. Name of the final annotated table file.
#'  Default is "anno.tsv".
#' @param tmpfold Character. Path to a temporary folder for intermediate files.
#'  Default is "./tmp".
#' @param saveobjectpath Character. Path to save intermediate objects. Default
#'  is NA.
#' @param savefinaltable Logical. Whether to save the final table to disk.
#'  Default is TRUE.
#' @param reload Logical. Whether to reload intermediate objects if available.
#'  Default is FALSE.
#' @param showtime Logical. Whether to display timing information. Default is
#'  FALSE.
#' @param showmemory Logical. Whether to display memory usage information.
#'  Default is FALSE.
#' @param deletetmp Logical. Whether to delete temporary files after processing.
#'  Default is TRUE.
#' @param verbose Logical. Whether to display detailed progress messages.
#'  Default is TRUE.
#'
#' @return A data frame representing the final table containing transcript
#' information and scores on 'windsize' windows for all experiments defined in
#' the experiment table (exptabpath).
#'
#' @details
#' The `preprocessing` function performs several key tasks:
#' 1. Filters Gencode annotations to retrieve "transcript" annotations.
#' 2. Differentiates between protein-coding (MANE_Select) and long non-coding
#'  (lncRNA, Ensembl_canonical) transcripts.
#' 3. Splits transcripts into windows of size `windsize`.
#' 4. Processes bedgraph files to retrieve values, exclude blacklisted regions,
#'  and retain high-mappability intervals.
#' 5. Generates a final annotated table with scores derived from the above
#' steps.
#'
#' Temporary files created during processing are optionally deleted at the end.
#'
#' @examples
#' # Example usage of preprocessing:
#' preprocessing(
#'   exptabpath = "./example_exptab.tsv",
#'   gencodepath = "./gencode.v38.annotation.gtf",
#'   windsize = 200,
#'   maptrackpath = "./mappability_track.bed",
#'   blacklistshpath = "./blacklist.bed",
#'   genomename = "hg38",
#'   nbcputrans = 2,
#'   finaltabpath = "./results/",
#'   finaltabname = "final_annotated_table.tsv",
#'   tmpfold = "./tmp",
#'   saveobjectpath = "./saved_objects",
#'   savefinaltable = TRUE,
#'   reload = FALSE,
#'   showtime = TRUE,
#'   showmemory = TRUE,
#'   deletetmp = TRUE,
#'   verbose = TRUE
#' )
#'
#' @seealso
#' [retrieveanno][makewindows][blacklisthighmap][createtablescores]
#' 
#' @export

preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans = 1, finaltabpath = "./",
    finaltabname = "anno.tsv", tmpfold = "./tmp", saveobjectpath = NA,
    savefinaltable = TRUE, reload = FALSE, showtime = FALSE, showmemory = FALSE,
    deletetmp = TRUE, verbose = TRUE) {

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
    .createbedgraphlistwmean(maptrackpath, blacklistshpath,
        exptabpath, nbcputrans, allwindowsbed, windsize, genomename, showtime,
        showmemory, saveobjectpath, reload, tmpfold, verbose)

    ## Creating the final table from the information retrieved from
    ## blacklisthighmap
    finaltable <- createtablescores(tmpfold, exptabpath, showmemory, showtime,
        savefinaltable, finaltabpath, finaltabname, verbose)

    ## Removing temporary files
    if (deletetmp) {
        if (verbose) message("\n ## Removing temporary files")
        unlink(tmpfold, recursive = TRUE)
    }

    if (showtime) {
        end_time_preprocessing <- Sys.time()
        timing <- end_time_preprocessing - start_time_preprocessing
            message("\n\n\t ## Total preprocessing in: ",
                format(timing, digits = 2))
    }

    return(finaltable)
}
