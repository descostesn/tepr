preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, nbcputrans = 1, nbcpubg = 1, finaltabpath = "./",
    finaltabname = "anno.tsv", saveobjectpath = NA, savefinaltable = TRUE,
    reload = FALSE, showstats = FALSE, showtime = FALSE, verbose = TRUE) {

    ## This function filters gencode annotations to retrieve "transcript". It
    ## then distinguishes transcripts coming from protein coding genes
    ## (MANE_Select) and those coming from long non-coding genes (lncRNA,
    ## Ensembl_canonical).
    if (verbose) message("## Filtering gencode annotations ##\n")
    if (!reload) {
        allannobed <- retrieveanno(exptabpath, gencodepath, saveobjectpath,
            showstats, verbose)
    } else {
        allannobedobj <- readRDS(file.path(saveobjectpath, "allannobed.rds"))
        if (verbose) message("Loading ", allannobedobj)
        allannobed <- readRDS(allannobedobj)
    }

    ## This functions uses the annotations filtered from gencode (allannobed).
    ## It removes any ensembl names containing "PAR_Y", filters out intervals
    ## smaller than windsize and splits each transcript into "windsize" windows.
    if (verbose) message("\n ## Splitting transcripts into windows ##\n")
    allwindowsbedobjpath <- file.path(saveobjectpath, "allwindowsbed.rds")
    if (!reload && !file.exists(allwindowsbedobjpath)) {
        allwindowsbed <- makewindows(allannobed, windsize, nbcputrans, verbose,
            saveobjectpath, showtime)
    } else {
        allwindowsbedobj <- readRDS(allwindowsbedobjpath)
        if (verbose) message("Loading ", allwindowsbedobj)
        allannobed <- readRDS(allwindowsbed)
    }

    ## Retrieving the values of the bedgraph files, removing black lists and
    ## keeping scores landing on high mappability intervals
    if (verbose) message("\n ## Retrieving the values of the bedgraph files, ",
        "removing black lists and keeping scores landing on high mappability",
        " intervals ##\n")
    bedgraphlistwmean <- blacklisthighmap(maptrackpath, blacklistshpath,
        exptabpath, nbcputrans, allwindowsbed, windsize, saveobjectpath, reload,
        showtime, verbose)

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

    return(finaltable)
}
