library("rtracklayer")
library("GenomicRanges")

##################
# PARAMETERS
##################

bgvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score"

allbgnicpath <- "/g/romebioinfo/tmp/preprocessing/backup/bedgraphwmeanlist.rds"
allwindowspath <- "/g/romebioinfo/tmp/preprocessing/allwindowsbed.rds"

blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/tmp/preprocessing/maptrackgr.gr"

##################
#FUNCTIONS
##################

bedtogr <- function(currentbed, strand = TRUE, symbol = TRUE,
    allwindows = FALSE) {

    if (!allwindows) {
        grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
                ranges = IRanges::IRanges(start = currentbed[, 2],
                                      end = currentbed[, 3],
                                      names = currentbed[, 4]),
                strand = if (strand) currentbed[, 6] else "*",
                symbol = if (symbol) currentbed[, 5] else NA)
    } else {
        grres <- GenomicRanges::GRanges(seqnames = currentbed[, 2],
                ranges = IRanges::IRanges(start = currentbed[, 3],
                                      end = currentbed[, 4],
                                      name = currentbed[, 5]),
                gene = currentbed[, 6], strand = currentbed[, 7],
                biotype = currentbed[, 1], window = currentbed[, 8],
                coord = currentbed[, 9])
    }
    return(grres)
}

##################
# MAIN
##################

## This is the ctrl rep1 fwd
bgvic <- read.delim(bgvicpath, header = FALSE)

## Selecting ctrl rep1 fwd
allbgnic <- readRDS(allbgnicpath)
names(allbgnic) <- gsub(".bg","",basename(names(allbgnic)))
bgnic <- allbgnic[["ctrl_rep1.forward"]]

## Reading all windows bed
allwindowsbed <- readRDS(allwindowspath)

## Selecting the lines corresponding to the gene ARF5
bgvicarf <- bgvic[which(bgvic$V6 == "ARF5"), ]
bgnicarf <- bgvic[which(bgnic$gene == "ARF5"), ]
allwindarf <- allwindowsbed[which(allwindowsbed$gene == "ARF5"), ]
allwindowsgr <- bedtogr(allwindarf, allwindows = TRUE)

## Reading black list and maptrack
blacklistbed <- read.delim(blacklistshpath, header = FALSE)
blacklistgr <- bedtogr(blacklistbed, strand = FALSE, symbol = FALSE)
maptrackgr <- readRDS(maptrackpath)
