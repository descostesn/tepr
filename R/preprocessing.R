## For different strings provided in the vector "valvec", perform the filtering
## of gentab on each string using the result of the previous filtering
.grepsequential <- function(valvec, gentab, invert = FALSE, verbose = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ] ## This line enables sequential grep
    }))
    return(gentab)
}

.sortedbedformat <- function(gencode) {
    gencode <- gencode[order(gencode$V1, gencode$V4), ] # nolint ## Ordering by chrom and start
    infolist <- strsplit(gencode$V9, ";")
    namevec <- gsub(" gene_name ", "", sapply(infolist, "[", 4)) # nolint
    ensnamevec <- gsub(" transcript_id ", "", sapply(infolist, "[", 2)) # nolint
    gencodebed <- cbind(gencode[, c(1, 4, 5)], ensnamevec, namevec,
        gencode[, 7])
    colnames(gencodebed) <- c("chrom", "start", "end", "ensembl", "symbol",
        "strand")
    return(gencodebed)
}


retrieveanno <- function(exptabpath, gencodepath, saveobjectpath = NA,
    verbose = TRUE) {

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Reading the information about experiments
    if (verbose) message("Reading the information about experiments")
    exptab <- read.csv(exptabpath, header = TRUE)
    checkexptab(exptab) # nolint

    ## Reading gencode file
    if (verbose) message("Reading gencode file and filtering")
    gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

    ## Keeping "transcript" annotations
    if (verbose) message("\t Keeping 'transcript' annotations")
    gencode <- gencode[which(gencode$V3 == "transcript"), ]

    ## Selecting Ensembl_canonical transcripts i.e. most representative
    ## transcript of the protein coding gene. This will be the MANE_Select
    ## transcript if there is one, or a transcript chosen by an Ensembl
    ## algorithm otherwise.
    if (verbose) message("\t Selecting Ensembl_canonical transcripts and ",
        "sorting")
    gencodeprotcod <- .grepsequential("MANE_Select", gencode)
    protcodbed <- .sortedbedformat(gencodeprotcod)

    ## Selecting long non-coding transcripts
    if (verbose) message("\t Selecting long non-coding transcripts and ",
        "sorting")
    lncrna <- .grepsequential(c("lncRNA", "Ensembl_canonical"), gencode)
    removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
    lncrna <- .grepsequential(removevec, lncrna, invert = TRUE)
    lncrnabed <- .sortedbedformat(lncrna)

}