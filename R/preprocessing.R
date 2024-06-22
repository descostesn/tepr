####################
# This script aims at performing pre-processing steps using R only.
#
# Descostes - June 2024 - R-4.4.1
####################


##################
# PARAMETERS
##################

gencodepath <- "/g/romebioinfo/Projects/tepr/downloads/gencode.v43.basic.annotation.gtf" # nolint


##################
#FUNCTIONS
##################

bedformat <- function(gencode) {
    gencode <- gencode[order(gencode$V1, gencode$V4), ] # nolint ## Ordering by chrom and start
    infolist <- strsplit(gencode$V9, ";")
    namevec <- gsub(" gene_name ", "", sapply(infolist, "[", 4)) # nolint
    ensnamevec <- gsub(" transcript_id ", "", sapply(infolist, "[", 2)) # nolint
    gencodebed <- cbind(gencode[, c(1, 4, 5)], ensnamevec, namevec,
        gencode[, 7])
    return(gencodebed)
}

##################
# MAIN
##################

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

## Selecting Ensembl_canonical transcripts i.e. most representative transcript
## of the gene. This will be the MANE_Select transcript if there is one, or a
## transcript chosen by an Ensembl algorithm otherwise.
gencode <- gencode[which(gencode$V3 == "transcript"), ]
gencode <- gencode[grep("MANE_Select", gencode$V9), ]

## Creating sorted bed format gencode
gencodebed <- bedformat(gencode)

