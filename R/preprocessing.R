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

.createbedgraphlistwmean <- function(maptrackpath, blacklistpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, genomename, showtime, showmemory,
    saveobjectpath, reload, tmpfold, chromtab, verbose) {

        if (verbose) message("\n ## Retrieving the values of the bedgraph ",
            "files, removing black lists and keeping scores landing on high ",
            "mappability intervals ##\n")
        blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
            nbcputrans, allwindowsbed, windsize, genomename, saveobjectpath,
            tmpfold, reload, showtime, showmemory, chromtab, forcechrom = TRUE,
            verbose)
}



#' Preprocess Experimental Data for Genomic Analysis
#'
#' @description
#' This function orchestrates a pipeline for preprocessing genomic data,
#' including filtering annotations, splitting transcripts into windows,
#' retrieving bedgraph values, and generating a final annotated table.
#'
#' @usage
#' preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
#' blacklistpath, genomename = NA, nbcputrans = 1, finaltabpath = tempdir(),
#' finaltabname = "anno.tsv", tmpfold = file.path(tempdir(), "tmptepr"),
#' saveobjectpath = tempdir(), savefinaltable = TRUE, reload = FALSE,
#' showtime = FALSE, showmemory = FALSE, deletetmp = TRUE, chromtab = NA,
#' forcechrom = FALSE, verbose = TRUE)
#'
#' @param exptabpath Character. Path to the experiment table file.
#' @param gencodepath Character. Path to the Gencode annotation file.
#' @param windsize Integer. Window size for splitting transcripts.
#' @param maptrackpath Character. Path to the mappability track file.
#' @param blacklistpath Character. Path to the blacklist file.
#' @param genomename Character. Name of the genome assembly (e.g., "hg38").
#'  Default is \code{NA}. If left to NA, chromtab should be provided.
#' @param nbcputrans Integer. Number of CPUs to use for transcript processing.
#'  Default is \code{1}.
#' @param finaltabpath Character. Path where the final annotated table will be
#'  saved. Default is \code{tempdir()}.
#' @param finaltabname Character. Name of the final annotated table file.
#'  Default is \code{anno.tsv}.
#' @param tmpfold Character. Path to a temporary folder for intermediate files.
#'  Default is \code{file.path(tempdir(), "tmptepr")}.
#' @param saveobjectpath Character. Path to save intermediate objects. Default
#'  is \code{tempdir()}.
#' @param savefinaltable Logical. Whether to save the final table to disk.
#'  Default is \code{TRUE}.
#' @param reload Logical. Whether to reload intermediate objects if available.
#'  Default is \code{FALSE}.
#' @param showtime Logical. Whether to display timing information. Default is
#'  \code{FALSE}.
#' @param showmemory Logical. Whether to display memory usage information.
#'  Default is \code{FALSE}.
#' @param deletetmp Logical. Whether to delete temporary files after processing.
#'  Default is \code{TRUE}.
#' @param chromtab A Seqinfo object retrieved with the rtracklayer method
#' \code{SeqinfoForUCSCGenome}. If NA, the method is called automatically and
#' the \code{genomename} should be provided. Default is \code{NA}.
#' @param forcechrom Logical indicating if the presence of non-canonical
#' chromosomes in chromtab (if not NA) should trigger an error. Default is
#' \code{FALSE}.
#' @param verbose Logical. Whether to display detailed progress messages.
#'  Default is \code{TRUE}.
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
#' \donttest{
#' ## Data
#' exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package = "tepr")
#' gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
#' maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
#'   package = "tepr")
#' blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
#'     package = "tepr")
#' windsize <- 200
#' genomename <- "hg38"
#'
#' ## Copying bedgraphs to the current directory
#' expdfpre <- read.csv(exptabpath)
#' bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
#'     package = "tepr"))
#' expdfpre$path <- bgpathvec
#' write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
#'     quote = FALSE)
#' exptabpath <- "exptab-preprocessing.csv"
#'
#' ## Testing preprocessing
#' finaltabtest <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
#'     blacklistpath, genomename = genomename, verbose = FALSE)}
#'
#' @seealso
#' [retrieveanno], [makewindows], [blacklisthighmap], [createtablescores]
#'
#' @export

preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = NA, nbcputrans = 1, finaltabpath = tempdir(),
    finaltabname = "anno.tsv", tmpfold = file.path(tempdir(), "tmptepr"),
    saveobjectpath = tempdir(), savefinaltable = TRUE, reload = FALSE,
    showtime = FALSE, showmemory = FALSE, deletetmp = TRUE, chromtab = NA,
    forcechrom = FALSE, verbose = TRUE) {
    
    if (!isTRUE(all.equal(typeof(chromtab), "S4"))) {
        if (is.na(genomename) && is.na(chromtab))
            stop("\n\t Either the genome name or chromtab should be ",
                "provided.\n")
    }

    if (reload && file.exists(file.path(saveobjectpath, "finaltable.rds")))
        stop("\n\t The final table already exists, set reload = FALSE to ",
            "create it again.\n")

    if (!is.na(chromtab) && !forcechrom) {
        if (!isTRUE(all.equal(is(chromtab), "Seqinfo")))
            stop("\n Chromtab should be a Seqinfo object. Use ",
                "rtracklayer::SeqinfoForUCSCGenome(genomename).\n")
        
        allchromvec <- GenomeInfoDb::seqnames(chromtab)
        idx <- grep("_|chrM", allchromvec, perl = TRUE, invert = FALSE)
        if (!isTRUE(all.equal(length(idx), 0)))
            stop("\n Non-canonical chromosomes found in chromtab. If you are ",
                "sure you want to proceed set forcechrom = TRUE.\n\n",
                paste(allchromvec[idx], collapse = " "))
    }

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
    .createbedgraphlistwmean(maptrackpath, blacklistpath,
        exptabpath, nbcputrans, allwindowsbed, windsize, genomename, showtime,
        showmemory, saveobjectpath, reload, tmpfold, chromtab, verbose)

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
