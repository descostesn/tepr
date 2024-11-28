preprocessing <- function(exptabpath, gencodepath, saveobjectpath = NA,
    verbose = TRUE) {

    ## This function filters gencode annotations to retrieve "transcript". It
    ## then distinguishes transcripts coming from protein coding genes
    ## (MANE_Select) and those coming from long non-coding genes (lncRNA,
    ## Ensembl_canonical).
    if (verbose) message("## Filtering gencode annotations\n")
    anno <- retrieveanno(exptabpath, gencodepath, saveobjectpath, verbose)
}