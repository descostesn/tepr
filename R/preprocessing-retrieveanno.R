## For different strings provided in the vector "valvec", perform the filtering
## of gentab on each string using the result of the previous filtering
.grepsequential <- function(valvec, gentab, verbose, invert = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message("\t\t", val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ] ## This line enables sequential grep
    }))
    return(gentab)
}

.sortedbedformat <- function(gencode) {
    ## Ordering by chrom and start
    gencode <- gencode[order(gencode$V1, gencode$V4), ] # nolint

    infolist <- strsplit(gencode$V9, ";")
    namevec <- gsub(" gene_name ", "", sapply(infolist, "[", 4)) # nolint
    ensnamevec <- gsub(" transcript_id ", "", sapply(infolist, "[", 2)) # nolint
    gencodebed <- cbind(gencode[, c(1, 4, 5)], ensnamevec, namevec,
        gencode[, 7])
    colnames(gencodebed) <- c("chrom", "start", "end", "ensembl", "symbol",
        "strand")
    return(gencodebed)
}


#' Retrieve and Combine Annotation Information
#'
#' @description
#' This function filters gencode annotations to retrieve "transcript". It then
#' distinguishes transcripts coming from protein coding genes (MANE_Select) and
#' those coming from long non-coding genes (lncRNA, Ensembl_canonical).
#'
#' @usage
#' retrieveanno(exptabpath, gencodepath, saveobjectpath = NA, showtime = FALSE,
#' verbose = TRUE)
#'
#' @param exptabpath Path to the experiment table file containing a table with
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param gencodepath Path to the GENCODE annotation file.
#' @param saveobjectpath Path to save intermediate R objects. Default is `NA`
#'  and R objects are not saved.
#' @param showtime Logical. If `TRUE`, displays timing information. Default is
#'  `FALSE`.
#' @param verbose Logical. If `TRUE`, provides detailed messages during
#'  execution. Default is `TRUE`.
#'
#' @return A data frame containing the combined annotation information for
#'  protein-coding and long non-coding RNA transcripts. If `saveobjectpath` is
#'  not `NA`, the object is also saved as an RDS file in the specified
#'  directory.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads experimental data from the provided CSV file and validates it.
#' 2. Reads genomic annotations from the gencode file and filters for
#'  transcripts.
#' 3. Separately processes protein-coding and long non-coding RNA transcripts:
#'    - For protein-coding genes, selects the most representative (MANE_Select
#'  or Ensembl_canonical) transcripts.
#'    - For long non-coding RNAs, filters out transcripts with undesirable
#'  evidence levels.
#' 4. Combines these annotations into a single data frame, labeling each
#'  transcript with its biotype.
#' 5. Optionally saves the resulting data frame as an RDS file in the specified
#'  directory.
#' 6. Optionally reports the total time taken for analysis.
#'
#' @examples
#' exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
#' gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
#'
#' ## Testing retrieveanno
#' allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
#'
#' @importFrom utils read.csv read.delim
#'
#' @export

retrieveanno <- function(exptabpath, gencodepath, saveobjectpath = NA,
    showtime = FALSE, verbose = TRUE) {

    if (showtime) start_time <- Sys.time()

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Reading the information about experiments
    if (verbose) message("Reading the information about experiments")
    exptab <- utils::read.csv(exptabpath, header = TRUE)
    checkexptab(exptab) # nolint

    ## Reading gencode file
    if (verbose) message("Reading gencode file and filtering")
    gencode <- utils::read.delim(gencodepath, header = FALSE, skip = 5)

    ## Keeping "transcript" annotations
    if (verbose) message("\t Keeping 'transcript' annotations")
    gencode <- gencode[which(gencode$V3 == "transcript"), ]

    ## Selecting Ensembl_canonical transcripts i.e. most representative
    ## transcript of the protein coding gene. This will be the MANE_Select
    ## transcript if there is one, or a transcript chosen by an Ensembl
    ## algorithm otherwise.
    if (verbose) message("\t Selecting Ensembl_canonical transcripts and ",
        "sorting")
    gencodeprotcod <- .grepsequential("MANE_Select", gencode, verbose)
    protcodbed <- .sortedbedformat(gencodeprotcod)

    ## Selecting long non-coding transcripts
    if (verbose) message("\t Selecting long non-coding transcripts and ",
        "sorting")
    lncrna <- .grepsequential(c("lncRNA", "Ensembl_canonical"), gencode,
        verbose)
    removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
    lncrna <- .grepsequential(removevec, lncrna, verbose, invert = TRUE)
    lncrnabed <- .sortedbedformat(lncrna)

    ## Combine the annotations
    if (verbose) message("\t Combine the annotations")
    protcodbed <- cbind(protcodbed, biotype = "protein-coding")
    lncrnabed <- cbind(lncrnabed, biotype = "lncRNA")
    allannobed <- rbind(protcodbed, lncrnabed)

    if (!is.na(saveobjectpath)) {
        outfile <- file.path(saveobjectpath, "allannobed.rds")
        if (verbose) message("\t Saving ", outfile)
        saveRDS(allannobed, outfile)
    }

    if (showtime) {
      end_time <- Sys.time()
      timing <- end_time - start_time
      message("\t\t ## Analysis performed in: ", format(timing, digits = 2))
    }

    return(allannobed)
}
