preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, nbcputrans = 1, nbcpubg = 1, saveobjectpath = NA,
    verbose = TRUE) {

    ## This function filters gencode annotations to retrieve "transcript". It
    ## then distinguishes transcripts coming from protein coding genes
    ## (MANE_Select) and those coming from long non-coding genes (lncRNA,
    ## Ensembl_canonical).
    if (verbose) message("## Filtering gencode annotations ##\n")
    allannobed <- retrieveanno(exptabpath, gencodepath, saveobjectpath, verbose)

    ## This functions uses the annotations filtered from gencode (allannobed).
    ## It removes any ensembl names containing "PAR_Y", filters out intervals
    ## smaller than windsize and splits each transcript into "windsize" windows.
    if (verbose) message("\n ## Splitting transcripts into windows ##\n")
    allwindowsbed <- makewindows(allannobed, windsize, nbcputrans, verbose,
        saveobjectpath)

    ## Retrieving the values of the bedgraph files, removing black lists and
    ## keeping scores landing on high mappability intervals
    if (verbose) message("\n ## Retrieving the values of the bedgraph files, ",
        "removing black lists and keeping scores landing on high mappability",
        " intervals ##\n")
    bedgraphlistwmean <- blacklisthighmap(maptrackpath, blacklistshpath,
        exptabpath, nbcputrans, allwindowsbed, windsize, saveobjectpath,
        verbose)

    ## Creating the final table from the information retrieved from
    ## blacklisthighmap
    if (verbose) message("\n ## Merging results of each bedgraph into a ",
        "single table ##\n")
    finaltable <- createtablescores(bedgraphlistwmean, nbcpubg, saveobjectpath,
        verbose)

    return(finaltable)
}