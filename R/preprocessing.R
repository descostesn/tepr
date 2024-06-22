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

grepsequential <- function(tokeepvec, gentab, invert = FALSE) {
    invisible(sapply(tokeepvec, function(tokeep) {
        idx <- grep(tokeep, gentab$V9)
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ]
    }))
    return(gentab)
}

sortedbedformat <- function(gencode) {
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
gencode <- gencode[which(gencode$V3 == "transcript"), ]

## Selecting Ensembl_canonical transcripts i.e. most representative transcript
## of the protein coding gene. This will be the MANE_Select transcript if there
## is one, or a transcript chosen by an Ensembl algorithm otherwise.
gencodeprotcod <- grepsequential("MANE_Select", gencode)
protcodbed <- sortedbedformat(gencodeprotcod)

## Retrieve long non-coding transcripts
lncrna <- grepsequential(c("lncRNA", "Ensembl_canonical"), gencode)
removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5", lncrna, invert = TRUE),
                grep("transcript_support_level 4", lncrna, invert = TRUE))



awk 'OFS="\t" {print $0}' Ensembl_canonical_TSL123.lncRNA.gtf | tr -d '";' | sort -k1,1 -k2,2n  | awk -F \t -v OFS='\t' '{print $0}' | awk -v OFS="\t" '{ print $1,$4,$5,$12,$16,$7}' > Ensembl_canonical_TSL123.lncRNA.bed



