library("rtracklayer")
library("GenomicRanges")
library("tibble")
library("valr")


##################
# PARAMETERS
##################

bgvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score"

allbgnicpath <- "/g/romebioinfo/tmp/preprocessing/backup/bedgraphwmeanlist.rds"
allwindowspath <- "/g/romebioinfo/tmp/preprocessing/allwindowsbed.rds"

exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv" # nolint
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint

nbcpubg <- 1



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

## Reading exptab, black list, and maptrack
exptab <- read.csv(exptabpath, header = TRUE)
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)

blacklistbed <- read.delim(blacklistshpath, header = FALSE)
# blacklistgr <- bedtogr(blacklistbed, strand = FALSE, symbol = FALSE)

maptrackbed <- read.delim(maptrackpath, header = FALSE)


## Debugging filtering

bedgraphgrlist <- retrieveandfilterfrombg(exptab, blacklistbed,
    maptrackbed, nbcpubg, expnamevec)



retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed, nbcpubg,
    allwindowsbed, expnamevec, verbose = TRUE) {

    if (verbose) message("Converting annotations' windows to tibble")
    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    ## Looping on each experiment bw file
    # currentpath <- exptab$path[1]
    # currentname <- expnamevec[1]
    # currentstrand <- exptab$strand[1]
    bedgraphgrlist <- parallel::mcmapply(function(currentpath, currentname,
        currentstrand, allwindowsbed, blacklistbed, maptrackbed, verbose) {

        if (verbose) message("\t Retrieving values for ", currentname)
        valgr <- rtracklayer::import.bedGraph(currentpath)
        if (verbose) message("\t\t Converting to tibble")
        valdf <- as.data.frame(valgr)
        colnames(valdf) <- c("chrom", "start", "end", "width", "strand",
            "score")
        valtib <- tibble::as_tibble(valdf)
        
        valr::bed_intersect(valtib, allwindtib)
!!!!!!!!!!!!!!!!


if (verbose) message("\t Keeping scores on annotations")
pairs <- IRanges::findOverlapPairs(valgr, allwindowsgr[which(strand(allwindowsgr) == "+")], ignore.strand = TRUE)
ansgr <- IRanges::pintersect(pairs, ignore.strand = TRUE)

if (verbose) message("\t\t Filtering out scores in black list ranges")
ansnoblackgr <- IRanges::subsetByOverlaps(ansgr, blacklistgr, invert = TRUE, ignore.strand = TRUE)

pairs <- IRanges::findOverlapPairs(ansnoblackgr, maptrackgr, ignore.strand = TRUE)
ansmaphigh <- IRanges::pintersect(pairs, ignore.strand = TRUE)


!!!!!!!!!!!!!!!!!



        if (verbose) message("\t\t Filtering out scores in black list ranges")
        resblack <- GenomicRanges::findOverlaps(valgr, blacklistgr,
            ignore.strand = TRUE)
        idxblack <- unique(S4Vectors::queryHits(resblack))
        BiocGenerics::score(valgr)[idxblack] <- NA

        if (verbose) message("\t\t Keeping high mappability scores")
        reshigh <- GenomicRanges::findOverlaps(valgr, maptrackgr,
            ignore.strand = TRUE)
        idxhigh <- unique(S4Vectors::queryHits(reshigh))
        if (isTRUE(all.equal(length(idxhigh), length(valgr))))
            message("Only highly mappable element were found")
        else
            ## Setting the scores of the ranges NOT in idxhigh to NA
            BiocGenerics::score(valgr)[-idxhigh] <- NA

        return(valgr)

    }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindowsbed,
        blacklistbed, maptrackbed, verbose), mc.cores = nbcpubg,
        SIMPLIFY = FALSE)

    return(bedgraphgrlist)
}



\