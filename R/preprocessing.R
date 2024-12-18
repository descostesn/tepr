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

# !!!!!!
# !!!!!! DOC TO DO
# !!!!!!

preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans = 1, nbcpubg = 1,
    finaltabpath = "./", finaltabname = "anno.tsv", tmpfold = "./tmp",
    saveobjectpath = NA, savefinaltable = TRUE, reload = FALSE,
    showtime = FALSE, showmemory = FALSE, verbose = TRUE) {

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
    !!
    # finaltable <- createtablescores(bedgraphlistwmean, nbcpubg, exptabpath,
    #     saveobjectpath, verbose)

    # if (savefinaltable) {
    #     outfile <- file.path(finaltabpath, finaltabname)
    #     if (verbose) message("\n ## Saving the final table to ", outfile)
    #     write.table(finaltable, file = outfile, sep = "\t", quote = FALSE,
    #         row.names = FALSE, col.names = FALSE)
    # }

    # !!!!!!! remove all saved obj if set to true (must be the default)

    if (showtime) {
        end_time_preprocessing <- Sys.time()
        timing <- end_time_preprocessing - start_time_preprocessing
            message("\n\n\t ## Total preprocessing in: ", timing) # nolint
    }

    #return(finaltable)
}
