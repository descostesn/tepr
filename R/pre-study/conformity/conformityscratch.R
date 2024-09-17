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
#allwindowsgr <- bedtogr(allwindarf, allwindows = TRUE)

## Reading exptab, black list, and maptrack
exptab <- read.csv(exptabpath, header = TRUE)
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)

blacklistbed <- read.delim(blacklistshpath, header = FALSE)
# blacklistgr <- bedtogr(blacklistbed, strand = FALSE, symbol = FALSE)

maptrackbed <- read.delim(maptrackpath, header = FALSE)


## Debugging filtering

bedgraphgrlist <- retrieveandfilterfrombg(exptab, blacklistbed,
    maptrackbed, nbcpubg, expnamevec)



retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed,
    nbcpubg, allwindowsbed, expnamevec, verbose = TRUE) {

    if (verbose) message("Converting annotations' windows to tibble")
    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    if (verbose) message("Converting blacklist to tibble")
    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

    if (verbose) message("Converting mappability track to tibble")
    colnames(maptrackbed) <- c("chrom", "start", "end", "id", "mapscore")
    maptracktib <- tibble::as_tibble(maptrackbed)

    ## Looping on each experiment bw file
    bedgraphlist <- parallel::mcmapply(function(currentpath, currentname,
        currentstrand, allwindtib, blacklisttib, maptracktib, nbcpuchrom,
        verbose) {

        ## Dealing with bedgraph values
        if (verbose) message("\t Retrieving values for ", currentname)
        valgr <- rtracklayer::import.bedGraph(currentpath)
        if (verbose) message("\t\t Converting to tibble")
        valdf <- as.data.frame(valgr)
        colnames(valdf) <- c("chrom", "start", "end", "width", "strand",
            "score")
        valtib <- tibble::as_tibble(valdf)

        ## Overlapping scores with anno on correct strand
        if (verbose) message("\t Retrieving scores on annotations of strand ",
            currentstrand)
        allwindstrand <- allwindtib %>% dplyr::filter(strand == currentstrand) # nolint
        resanno <- valr::bed_intersect(valtib, allwindstrand,
            suffix = c("", ".anno"))

        ## Removing black list
        if (verbose) message("\t Keeping scores not on black list")
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)

        ## Keeping scores on high mappability track
        if (verbose) message("\t Keeping scores on high mappability track")
        chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
        resmaplist <- lapply(chromvec, function(currentchrom) {
            if (verbose) message("\t\t\t over ", currentchrom)
            resmap <-  tryCatch({
                valr::bed_intersect(resblack,
                maptracktib  %>% dplyr::filter(chrom == currentchrom), # nolint
                suffix = c("", "maphigh"))
            }, error = function(e) {
                message(e)
                return(NA)
            })
            return(resmap)})
        if (verbose) message("\t\t Merging chromosomes")
        resmap <- purrr::map_dfr(resmaplist, dplyr::bind_rows)

        return(res)

    }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, verbose), SIMPLIFY = FALSE,
        mc.cores = nbcpubg)

    return(bedgraphlist)
}











!!!!!!!!!!!!!!!

bedgraphwmeanreplace <- function(bedgraphgrlist, exptab, expnamevec,
    allwindowsgr, windsize, nbcputrans) {

        bedgraphwmeanlist <- mapply(function(currentgr, currentstrand,
            currentname, allwindowsgr, windsize, nbcputrans) {

                message("Overlapping ", currentname, " with annotations on ",
                    "strand ", currentstrand)
                BiocGenerics::strand(currentgr) <- currentstrand
                res <- GenomicRanges::findOverlaps(currentgr, allwindowsgr,
                    ignore.strand = FALSE)

                message("\t Building scoring results by transcript")
                ## Separating the bedgraph score indexes by transcript names
                idxanno <- S4Vectors::subjectHits(res)
                idxbgscorebytrans <- split(as.data.frame(res),
                    factor(names(allwindowsgr)[idxanno]))

                ## For each transcript, retrieve the information and the
                ## bedgraph coordinates, strand and scores, applying a weighted
                ## mean
                message("\t Weighted mean on duplicated frames for each ",
                    "transcript")
                start_time <- Sys.time()
                dfwmeanbytranslist <- summarizebywmean(idxbgscorebytrans,
                    allwindowsgr, currentgr, currentstrand, currentname,
                    windsize, nbcputrans)
                end_time <- Sys.time()
                message("\t\t ## Analysis performed in: ",
                    end_time - start_time)

                if (!isTRUE(all.equal(unique(sapply(dfwmeanbytranslist,nrow)),
                    windsize)))
                    stop("Problem in replacing scores by weighted mean",
                        " on the data")

                message("\t Combining transcripts")
                dfwmeanbytrans <- do.call("rbind", dfwmeanbytranslist)
                end_time2 <- Sys.time()
                message("\t\t ## Analysis performed in: ",
                    end_time2 - start_time)

                return(dfwmeanbytrans)

        }, bedgraphgrlist, exptab$strand, expnamevec,
            MoreArgs = list(allwindowsgr, windsize, nbcputrans),
            SIMPLIFY = FALSE)
        return(bedgraphwmeanlist)
}

