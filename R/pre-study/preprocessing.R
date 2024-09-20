####################
# This script aims at performing the pre-processing steps using R only.
#
# Descostes - June 2024 - R-4.4.1
####################


library("AnnotationHub")
library("GenomeInfoDb")
library("GenomicRanges")
library("rtracklayer")
library("parallel")
library("purrr")
library("dplyr")



##################
# PARAMETERS
##################

gencodepath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
windsize <- 200
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv" # nolint
database_name <- "org.Hs.eg.db"
outputfolder <- "/g/romebioinfo/tmp/preprocessing"
robjoutputfold <- outputfolder

## Set this variable to NULL if the online retrieval should be performed
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
## The bed file below was created and sent by Victor
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
## Parallelization on bedgraph files. The maximum should be equal to the number of bedgraph files.  # nolint
nbcpubg <- 8
## Parallelization on transcripts. The maximum should be limited to the capacity of your machine.  # nolint
nbcputrans <- 20

## Note: For a complete list of blacklist names see
## ah <- AnnotationHub() # nolint
## query_data <- subset(ah, preparerclass == "excluderanges") # nolint
## print(query_data) # nolint
blacklistname <- "hg38.Kundaje.GRCh38_unified_Excludable"



##################
#FUNCTIONS
##################


createfolder <- function(outfold) {
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
}

checkexptab <- function(exptab) {
    colnamevec <- c("condition", "direction", "path", "replicate", "strand")
    if (!isTRUE(all.equal(sort(colnames(exptab)), colnamevec)))
        stop("The experiment table should have the columns: ",
            "'condition', 'direction', 'path', 'replicate', 'strand'")

    if (!isTRUE(all.equal(length(unique(exptab$condition)), 2)))
        stop("The table should only contain two conditions")

    if (isTRUE(all.equal(length(grep("ctrl", exptab$condition)), 0)))
        stop("The control condition (or condition 1) should be designated",
            " by 'ctrl'")

    directionvec <- unique(exptab$direction)
    if (!isTRUE(all.equal(length(directionvec), 2)) ||
        !isTRUE(all.equal(directionvec, c("fwd", "rev"))))
        stop("Only two values are allowed for the column direction of the",
            "experiment table, 'fwd' and 'rev'")

    strandvec <- unique(exptab$strand)
    if (!isTRUE(all.equal(strandvec, c("+", "-"))))
        stop("The strand column of the experiment table should only contain",
            " '+' and '-'.")
}

createblacklist <- function(blacklistname, outputfolder) { # nolint

    blacklistgr <- AnnotationHub::query(AnnotationHub::AnnotationHub(),
        blacklistname)[[1]]
    blacklistgr <- blacklistgr %>%
                   sort() %>%
                   GenomeInfoDb::keepStandardChromosomes(pruning.mode = "tidy")
    createfolder(outputfolder)
    write.table(as.data.frame(blacklistgr),
            file = file.path(outputfolder, paste0(blacklistname, ".bed")),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(blacklistgr)
}

## For different strings provided in the vector "valvec", perform the filtering
## of gentab on each string using the result of the previous filtering
grepsequential <- function(valvec, gentab, invert = FALSE, verbose = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ] ## This line enables sequential grep
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
    colnames(gencodebed) <- c("chrom", "start", "end", "ensembl", "symbol",
        "strand")
    return(gencodebed)
}

checkremoval <- function(datagr, dataremovedgr, dataname, removename,
    toremovegr, removeopt) {

    message("\n\n Checking operations for ", dataname, " with data of ",
        removename)

    ## Calculating number of elements overlapping toremovegr
    res <- GenomicRanges::findOverlaps(datagr, toremovegr)
    nboverdata <- length(unique(S4Vectors::queryHits(res)))

    if (removeopt) {
        message("The number of elements that should be removed is: ",
            nboverdata)
        subres <- length(datagr) - nboverdata
        message("The number of elements of the resulting object after ",
            "subtraction should be: ", length(datagr), "-", nboverdata, "=",
            subres)
        message("The number of elements in the resulting object is: ",
            length(dataremovedgr))
    } else {
        message("The number of elements of the data to keep is: ",
            nboverdata)
        message("The number of elements before the overlap is: ",
            length(datagr))
        message("The number of elements in the resulting object is: ",
            length(dataremovedgr))
    }
}


.divideannoinwindows <- function(expbed, windcoordvec, nbwindows, nbcputrans) {

    cl <- parallel::makeCluster(nbcputrans)
    windflist <- parallel::parLapply(cl, seq_len(nrow(expbed)),
        function(i, expbed, windcoordvec, nbwindows) {

            currentanno <- expbed[i, ]
            ## Retrieve the necessary gene information
            currentstart <- currentanno$start
            currentend <- currentanno$end
            currentstrand <- currentanno$strand
            windowvec <- windcoordvec

            ## Compute the vector with the size of each window
            lgene <- currentend - currentstart
            windowsize <- round(lgene / nbwindows)
            missingbp <- lgene %% nbwindows
            windsizevec <- rep(windowsize, nbwindows)
            ## Add the missing nb of bp (that is ignore by tile) in the last
            ## element of windsizevec
            if (!isTRUE(all.equal(missingbp, 0)))
                windsizevec[nbwindows] <- windsizevec[nbwindows] + missingbp

            ## Building the start and end vectors using the cummulative sum
            cumsumvec <- cumsum(c(currentstart, windsizevec))
            startvec <- cumsumvec[-length(cumsumvec)]
            endvec <- cumsumvec[-1]
            if (!isTRUE(all.equal(endvec - startvec, windsizevec)))
                stop("Problem in the calculation of windows")

            ## Inverting start, end, and window vectors if strand is negative
            if (isTRUE(all.equal(currentstrand, "-"))) {
                startvec <- rev(startvec)
                endvec <- rev(endvec)
                windowvec <- rev(windcoordvec)
            }

            ## Build the result data.frame containing the coordinates of each
            ## frame alongside window and coord numbers
            res <- data.frame(biotype = currentanno$biotype,
                chr = currentanno$chrom, coor1 = startvec,
                coor2 = endvec,  transcript = currentanno$ensembl,
                gene = currentanno$symbol, strand = currentstrand,
                window = windowvec, coord = windcoordvec)
            return(res)}, expbed, windcoordvec, nbwindows)
    stopCluster(cl)
    nbwindcheck <- unique(sapply(windflist, nrow))
    if (!isTRUE(all.equal(length(nbwindcheck), 1)) ||
        !isTRUE(all.equal(nbwindcheck, 200)))
        stop("Problem in the nb of windows per transcript retrieved")
    windf <- do.call("rbind", windflist)

    return(windf)
}

makewindowsbedtools <- function(expbed, nbwindows, nbcputrans, verbose = TRUE) {

    ## Filtering out intervals smaller than nbwindows
    idxsmall <- which((expbed$end - expbed$start) < nbwindows)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        message("Excluding ", lsmall, "/", nrow(expbed), " annotations that ",
        "are too short.")
        expbed <- expbed[-idxsmall,]
    }

    ## Splitting each transcript into "nbwindows" windows
    if (verbose) message("\t Splitting ", nrow(expbed), " transcript into ",
        nbwindows, " windows data.frame")
    windcoordvec <- seq_len(nbwindows)
    winddf <- .divideannoinwindows(expbed, windcoordvec, nbwindows,
        nbcputrans)

    return(winddf)
}


.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed) {

    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

    colnames(maptrackbed) <- c("chrom", "start", "end", "id", "mapscore")
    maptracktib <- tibble::as_tibble(maptrackbed)

    return(list(allwindtib, blacklisttib, maptracktib))
}

.retrievebgval <- function(currentpath, verbose) {

    valgr <- rtracklayer::import.bedGraph(currentpath)
    if (verbose) message("\t\t Converting to tibble")
    valdf <- as.data.frame(valgr)
    colnames(valdf) <- c("chrom", "start", "end", "width", "strand", "score")
    valtib <- tibble::as_tibble(valdf)
    return(valtib)
}

.removeblacklist <- function(allwindstrand, valtib, currentstrand,
    blacklisttib) {
        ## Retrieving scores on annotations of strand
        suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
                suffix = c("", ".window")))
        ## Removing black list
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)
        return(resblack)
}

.retrieveonhighmap <- function(resblack, maptracktib, currentchrom) {

    ## Keeping scores on high mappability track
    resmap <-  valr::bed_intersect(resblack,
        maptracktib  %>% dplyr::filter(chrom == currentchrom), # nolint
        suffix = c(".bg", ".maphigh"))
    colnames(resmap) <- gsub(".window.bg", ".window", colnames(resmap))

    ## Removing mapping columns and duplicates
    resmap <- resmap[, -grep(".maphigh|.overlap|.source", colnames(resmap))]
    resmap <- resmap %>% dplyr::distinct(chrom, start.bg, end.bg, # nolint
        start.window, end.window, .keep_all = TRUE) # nolint
    invisible(gc())
    return(resmap)
}

.uniqueformatcolnames <- function(currenttrans) {

    ## Verifying uniformity of chrom, transcript, and genes
    uniquechrom <- as.character(unique(currenttrans$chrom))
    uniquetrans <- as.character(unique(currenttrans$transcript.window))
    uniquegene <- as.character(unique(currenttrans$gene.window))

    if (!isTRUE(all.equal(length(uniquechrom), 1)) ||
        !isTRUE(all.equal(length(uniquetrans), 1)) ||
        !isTRUE(all.equal(length(uniquegene), 1)))
            stop("chrom, transcript, and genes should be unique, this should", # nolint
                " not happen. Contact the developper.") # nolint

    ## Renaming window and coord columns removing the suffix
    colnamevec <- colnames(currenttrans)
    colnames(currenttrans)[which(colnamevec == "window.window")] <- "window"
    colnames(currenttrans)[which(colnamevec == "coord.window")] <- "coord"
    colnames(currenttrans)[which(colnamevec == "gene.window")] <- "gene"
    colnames(currenttrans)[which(colnamevec == "transcript.window")] <- "transcript" # nolint

    return(list(currenttrans, uniquechrom, uniquetrans, uniquegene))

}

.retrievemissingwind <- function(idxnavec, allwindstrand, currenttrans,
    uniquechrom, uniquetrans, uniquegene) {

    ## For each missing window whose number is contained in idxnavec
    missingrowslist <- lapply(idxnavec, function(idxna, allwindstrand) {

        ## Retrieving the line of the missing window in allwindstrand
        idxmissing <-  which(allwindstrand$chrom == uniquechrom &
            allwindstrand$transcript == uniquetrans &
            allwindstrand$gene == uniquegene &
            allwindstrand$window == idxna)

        if (!isTRUE(all.equal(length(idxmissing), 1)))
            stop("Problem in retrieving the missing window, this should not ",
                "happen. Contact the developper.")

        ## Below the bedgraph information columns are set to NA. These columns will be removed later # nolint
        ## The score is set to NA since it is a missing value resulting from removing black list and low mappability (keeping high mappability) # nolint
        ## Filling the other columns with the line retrieved in allwindstrand # nolint
        windstrandrow <- allwindstrand[idxmissing, ]
        resmissing <- data.frame(chrom = windstrandrow$chrom,
            start.bg = NA, end.bg = NA, width.bg = NA, strand.bg = "*", ## Set the bedgraph info # nolint
            score.bg = NA, ## Set the score to NA to keep track of missing values # nolint
            biotype.window = windstrandrow$biotype,
            start.window = windstrandrow$start,
            end.window = windstrandrow$end,
            transcript = windstrandrow$transcript, gene = windstrandrow$gene,
            strand.window = windstrandrow$strand, window = windstrandrow$window,
            coord = windstrandrow$coord)

        return(resmissing)
        }, allwindstrand)

    missingrowsdf <- do.call("rbind", missingrowslist)
    currenttrans <- rbind(currenttrans, missingrowsdf)
    currenttrans <- currenttrans[order(currenttrans$coord), ]

    return(currenttrans)
}

.arrangewindows <- function(currenttrans, windsize, allwindstrand,
    currentname) {

    res <- .uniqueformatcolnames(currenttrans)
    currenttrans <- res[[1]]
    uniquechrom <- res[[2]]
    uniquetrans <- res[[3]]
    uniquegene <- res[[4]]

    ## Identifying the missing window in currenttrans
    idx <- match(seq_len(windsize), unique(currenttrans$window))
    idxnavec <- which(is.na(idx))

    ## If some windows are missing
    if (!isTRUE(all.equal(length(idxnavec), 0)))
        currenttrans <- .retrievemissingwind(idxnavec,
            allwindstrand, currenttrans, uniquechrom,
            uniquetrans, uniquegene)

    idxscore <- which(colnames(currenttrans) == "score.bg")
    scorename <- paste0(currentname, "_score") # nolint
    colnames(currenttrans)[idxscore] <- scorename

    return(list(currenttrans, uniquetrans))
}

.computewmeanvec <- function(dupframenbvec, currenttrans, currentname,
    colscore) {

    ## For each duplicated frame
    wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans, currentname,
        colscore) {

        ## Selecting all rows having a window equal tro nbdup
        allframedf <- currenttrans[which(currenttrans$window == nbdup), ]
        if (isTRUE(all.equal(nrow(allframedf), 1)))
            stop("There should be more than one frame selected")

        ## Testing that the coord of the window is the same for all scores
        ## selected (this should not give an error)
        windowstart <- unique(allframedf$start.window)
        windowend <- unique(allframedf$end.window)
        if (!isTRUE(all.equal(length(windowstart), 1)) ||
            !isTRUE(all.equal(length(windowend), 1)))
                stop("The size of the window is not unique for the frame rows ",
                    "selected, this should not happen, contact the developper.")

        ## Retrieve the nb of overlapping nt for each score
        overntvec <- apply(allframedf, 1,
            function(x, currentname, windowstart, windowend) {
                nt <- seq(from = x["start.bg"], to = x["end.bg"], by = 1)
                overnt <- length(which(nt >= windowstart & nt <= windowend))
                return(overnt)
            }, currentname, windowstart, windowend)

        ## Computing weighted mean
        allscores <- as.data.frame(allframedf[,colscore])[[1]]
        wmean <- weighted.mean(allscores, overntvec)
        return(wmean)
    }, currenttrans, currentname, colscore)
    return(wmeanvec)
}

.replaceframeswithwmean <- function(currenttrans, dupidx, windsize, transname,
    dupframenbvec, colscore, wmeanvec) {

        ## Remove duplicated frames and replace scores by wmean
        currenttrans <- currenttrans[-dupidx, ]
        if (!isTRUE(all.equal(nrow(currenttrans), windsize)))
            stop("The number of frames should be equal to windsize: ",
                windsize, " for transcript ", transname)
        idxscorereplace <- match(dupframenbvec, currenttrans$window)
        if (!isTRUE(all.equal(dupframenbvec,
            currenttrans$window[idxscorereplace])))
            stop("Problem in replacing scores by wmean, contact the developer.")
        currenttrans[idxscorereplace, colscore] <- wmeanvec

        return(currenttrans)
}

.missingandwmean <- function(resmap, windsize, allwindstrand, currentname,
    nbcputrans) {

    ## Splitting the scores kept on the high mappability track by transcript
    bgscorebytrans <- split(resmap, factor(resmap$transcript.window))

    ## Performing identification of missing windows and calculation of weighted
    ## means for each transcript. This is parallelized on nbcputrans CPUs.
    bytranslist <- parallel::mclapply(bgscorebytrans, function(currenttrans,
        windsize, allwindstrand, currentname) {

            ## Identification of missing windows for the current transcript
            ## and setting their scores to NA. Indeed if a window is missing
            ## it is because it was in a black list or a low mappability
            ## interval.
            res <- .arrangewindows(currenttrans, windsize, allwindstrand,
                currentname)
            currenttrans <- res[[1]]
            transname <- res[[2]]

            ## Identifying duplicated windows that will be used to compute
            ## a weighted mean.
            dupidx <- which(duplicated(currenttrans$window))
            colscore <- paste0(currentname, "_score") # nolint

            if (!isTRUE(all.equal(length(dupidx), 0))) {
                dupframenbvec <- unique(currenttrans$window[dupidx])
                ## For each duplicated frame, compute the weighted mean
                wmeanvec <- .computewmeanvec(dupframenbvec, currenttrans,
                    currentname, colscore)
                ## Remove duplicated frames and replace scores by wmean
                currenttrans <- .replaceframeswithwmean(currenttrans, dupidx,
                    windsize, transname, dupframenbvec, colscore, wmeanvec)
             }

             currenttrans <- currenttrans[order(currenttrans$coord), ]
             return(currenttrans)
    }, windsize, allwindstrand, currentname, mc.cores = nbcputrans)

    transdf <- do.call("rbind", bytranslist)
    rm(bytranslist)
    invisible(gc())
    return(transdf)
}

.processingbychrom <- function(maptracktib, allwindstrand, currentname,
    resblack, nbcputrans, verbose) {

        if (verbose) message("\t\t Retrieving list of chromosomes")
        chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
        if (verbose) message("\t\t Formatting scores")
        bychromlist <- lapply(chromvec, function(currentchrom, allwindstrand,
            currentname, resblack, maptracktib, nbcputrans) {

                if (verbose) message("\t\t\t over ", currentchrom)

                ## Retrieving scores on high mappability intervals
                if (verbose) message("\t\t\t Keeping scores on high ",
                    "mappability track")
                resmap <- .retrieveonhighmap(resblack, maptracktib,
                    currentchrom)

                ## Processing data per transcript for windows and wmean
                message("\t\t\t Setting missing windows scores to NA and",
                    " computing weighted mean for each transcript")
                transdf <- .missingandwmean(resmap, windsize, allwindstrand,
                    currentname, nbcputrans)
                return(transdf)
            }, allwindstrand, currentname, resblack, maptracktib, nbcputrans)

            ## Merging results that were computed on each chromosome
            if (verbose) message("\t\t Merging results that were computed on",
                " each chromome")
            resallchrom <- do.call("rbind", bychromlist)
            rm(bychromlist)
            invisible(gc())
            return(resallchrom)
}

retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed, # nolint
    nbcputrans, allwindowsbed, expnamevec, windsize, verbose = TRUE) {

        if (verbose) message("\t Converting annotations' windows, blacklist,",
            " and mappability track to tibble")
        tibres <- .convertotibble(allwindowsbed, blacklistbed, maptrackbed)
        allwindtib <- tibres[[1]]
        blacklisttib <- tibres[[2]]
        maptracktib <- tibres[[3]]

        ## Looping on each experiment bg file
        if (verbose) message("\t For each bedgraph file")
        bedgraphlistwmean <- mapply(function(currentpath, currentname,
            currentstrand, allwindtib, blacklisttib, maptracktib, windsize,
            nbcputrans, verbose) {

            ## Retrieving bedgraph values
            if (verbose) message("\t\t Retrieving begraph values for ",
                currentname)
            valtib <- .retrievebgval(currentpath, verbose)

            ## Keeping window coordinates on the correct strand
            if (verbose) message("\t\t Retrieving coordinates on strand ",
                currentstrand)
            allwindstrand <- allwindtib %>%
                dplyr::filter(strand == as.character(currentstrand)) # nolint

            ## Overlapping scores with anno on correct strand and remove
            ## blacklist
            if (verbose) message("\t\t Keeping scores outside blacklist ",
                "intervals")
            resblack <- .removeblacklist(allwindstrand, valtib, currentstrand,
                blacklisttib)
            rm(valtib)
            invisible(gc())

            ## Processing by chromosomes because of size limits, the mappability
            ## track has too many rows. Formatting scores, keeping those on
            ## high mappability, filling missing windows, and compute wmean
            resallchrom <- .processingbychrom(maptracktib, allwindstrand,
                currentname, resblack, nbcputrans, verbose)
            return(resallchrom)
        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, nbcputrans, verbose),
        SIMPLIFY = FALSE)

        return(bedgraphlistwmean)
}


createtablescores <- function(bedgraphlistwmean, nbcpubg) {

    ## Creating a rowid that will be used for merging
    message("\t Adding rowid for each bedgraph") # nolint
    rowidreslist <- parallel::mclapply(bedgraphlistwmean, function(tab) {

        rowidvec <- paste(tab$biotype.window, tab$chrom, tab$start.window,
            tab$end.window, tab$strand.window, tab$gene,
            tab$transcript, paste0("frame", tab$window),
            paste0("coord", tab$coord), sep = "_")

        ## Inserting rowid col after transcript
        tab <- tab %>% tibble::add_column(rowid = rowidvec,
            .after = "transcript")

        ## Remove bedgraph columns
        tab <- tab[, -grep(".bg", colnames(tab))]

        ## Move the score column at the end of the table
        colnamevec <- colnames(tab)
        idxscore <- grep("_score", colnamevec)
        if (!isTRUE(all.equal(length(idxscore), 1)))
            stop("When creating the final table, the score column is not ",
                "unique for a given bedgraph. This should not happen. Contact",
                " the developer.")
        tab <- dplyr::relocate(tab, colnamevec[idxscore], .after = "coord")
        return(tab)
    }, mc.cores = nbcpubg)

    message("\t Joining the elements of each bedgraph")
    completeframedf <- purrr::reduce(rowidreslist, dplyr::full_join,
        by = c("chrom", "start.window", "end.window", "strand.window", "gene",
        "biotype.window", "window", "coord", "transcript", "rowid"))

    return(completeframedf)
}




##################
# MAIN
##################

createfolder(robjoutputfold)

## Reading the information about experiments
exptab <- read.csv(exptabpath, header = TRUE)
checkexptab(exptab)

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

## Keeping "transcript" annotations
gencode <- gencode[which(gencode$V3 == "transcript"), ]

## Selecting Ensembl_canonical transcripts i.e. most representative transcript
## of the protein coding gene. This will be the MANE_Select transcript if there
## is one, or a transcript chosen by an Ensembl algorithm otherwise.
gencodeprotcod <- grepsequential("MANE_Select", gencode)
protcodbed <- sortedbedformat(gencodeprotcod)

## Retrieve long non-coding transcripts
lncrna <- grepsequential(c("lncRNA", "Ensembl_canonical"), gencode)
removevec <- c("not_best_in_genome_evidence", "transcript_support_level 5",
                "transcript_support_level 4")
lncrna <- grepsequential(removevec, lncrna, invert = TRUE)
lncrnabed <- sortedbedformat(lncrna)

## Combine the annotations
message("Combine the annotations")
protcodbed <- cbind(protcodbed, biotype = "protein-coding")
lncrnabed <- cbind(lncrnabed, biotype = "lncRNA")
allannobed <- rbind(protcodbed, lncrnabed)

saveRDS(allannobed, file.path(robjoutputfold, "allannobed.rds"))

## Make windows for all annotations
message("Make windows for all annotations")
idxpar <- grep("PAR_Y", allannobed$ensembl)
if (!isTRUE(all.equal(length(idxpar), 0)))
    allannobed <- allannobed[-idxpar, ]

allwindowsbed <- makewindowsbedtools(expbed = allannobed, nbwindows = windsize,
    nbcputrans = nbcputrans)
saveRDS(allwindowsbed, file.path(robjoutputfold, "allwindowsbed.rds"))


## Retrieving the values of the bedgraph files, removing black lists and keeping
## high mappability scores
message("Reading the black list and mappability track")
if (is.null(blacklistshpath)) {
    message("Retrieving the black list online")
    blacklistgr <- createblacklist(blacklistname, outputfolder)
} else {
    blacklistbed <- read.delim(blacklistshpath, header = FALSE)
}

maptrackbed <- read.delim(maptrackpath, header = FALSE)
saveRDS(maptrackbed, file.path(robjoutputfold, "maptrackbed.rds"))

message("Removing scores within black list intervals, keeping those on high",
" mappability regions, and computing weighted means.")
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
bedgraphlistwmean <- retrieveandfilterfrombg(exptab, blacklistbed, maptrackbed,
    nbcputrans, allwindowsbed, expnamevec, windsize)

message("Merging results of each bedgraph into a single table")
finaltab <- createtablescores(bedgraphlistwmean, nbcpubg)
saveRDS(finaltab, file = file.path(outputfolder, "finaltab.rds"))
