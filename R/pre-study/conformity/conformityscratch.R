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
windsize <- 200

##################
#FUNCTIONS
##################



##################
# MAIN
##################

## Reading all windows bed
allwindowsbed <- readRDS(allwindowspath)

## Reading exptab, black list, and maptrack
exptab <- read.csv(exptabpath, header = TRUE)
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
blacklistbed <- read.delim(blacklistshpath, header = FALSE)
maptrackbed <- read.delim(maptrackpath, header = FALSE)


##################### TEST
## This is the ctrl rep1 fwd
bgvic <- read.delim(bgvicpath, header = FALSE)

## Selecting ctrl rep1 fwd
allbgnic <- readRDS(allbgnicpath)
names(allbgnic) <- gsub(".bg","",basename(names(allbgnic)))
bgnic <- allbgnic[["ctrl_rep1.forward"]]
rm(allbgnic)
gc()
## Selecting the lines corresponding to the gene ARF5
bgvicarf <- bgvic[which(bgvic$V6 == "ARF5"), ]
bgnicarf <- bgvic[which(bgnic$gene == "ARF5"), ]
allwindarf <- allwindowsbed[which(allwindowsbed$gene == "ARF5"), ]
#allwindowsgr <- bedtogr(allwindarf, allwindows = TRUE)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

.convertotibble <- function(allwindowsbed, blacklistbed, maptrackbed, verbose) {

    if (verbose) message("Converting annotations' windows to tibble") # nolint
    colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window", "coord")
    allwindtib <- tibble::as_tibble(allwindowsbed)

    if (verbose) message("Converting blacklist to tibble") # nolint
    colnames(blacklistbed) <- c("chrom", "start", "end", "type")
    blacklisttib <- tibble::as_tibble(blacklistbed)

    if (verbose) message("Converting mappability track to tibble") # nolint
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
    blacklisttib, verbose) {

        if (verbose) message("\t Retrieving scores on annotations of strand ", # nolint
                currentstrand)
        suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
                suffix = c("", ".window")))

        ## Removing black list
        if (verbose) message("\t Keeping scores not on black list") # nolint
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

retrieveandfilterfrombg <- function(exptab, blacklistbed, maptrackbed,
    nbcpubg, allwindowsbed, expnamevec, windsize, verbose = TRUE) {

        tibres <- .convertotibble(allwindowsbed, blacklistbed, maptrackbed,
            verbose)
        allwindtib <- tibres[[1]]
        blacklisttib <- tibres[[2]]
        maptracktib <- tibres[[3]]

        ## Looping on each experiment bg file
        bedgraphlist <- parallel::mcmapply(function(currentpath, currentname,
            currentstrand, allwindtib, blacklisttib, maptracktib, nbcpuchrom,
            windsize, verbose) {

            ## Retrieving bedgraph values
            if (verbose) message("\t Retrieving values for ", currentname) # nolint
            valtib <- .retrievebgval(currentpath, verbose)

            ## Keeping window coordinates on the correct strand
            allwindstrand <- allwindtib %>%
                dplyr::filter(strand == as.character(currentstrand)) # nolint

            ## Overlapping scores with anno on correct strand and remove
            ## blacklist
            resblack <- .removeblacklist(allwindstrand, valtib, currentstrand,
                blacklisttib, verbose)

            ## Processing by chromosomes because of size limits, the mappability
            ## track has too many rows
            if (verbose) message("\t Keeping scores on high mappability track")
            chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
            resmaplist <- lapply(chromvec, function(currentchrom, allwindstrand,
                currentname, resblack, maptracktib) {

                if (verbose) message("\t\t over ", currentchrom)

                if (verbose) message("\t\t\t Keeping scores on high ",
                    "mappability track")
                resmap <- .retrieveonhighmap(resblack, maptracktib,
                    currentchrom)

                ## Processing data per transcript
                message("\t\t\t Building scoring results by transcript")
                bgscorebytrans <- split(resmap,
                    factor(resmap$transcript.window))

                #currenttrans=bgscorebytrans[[1]]
                bytranslist <- lapply(bgscorebytrans,
                    function(currenttrans, windsize, allwindstrand,
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
                    if (!isTRUE(all.equal(length(idxnavec), 0))) {

                        .retrievemissingwind <- function(idxnavec, allwindstrand, currentstrand, currenttrans) {}

                        ## For each missing window whose number is contained in idxnavec
                        #message("\t\t\t Retrieving missing scores")
                        missingrowslist <- lapply(idxnavec, function(idxna, allwindstrand) {
                            ## Retrieving the line of the missing window in allwindstrand
                            idxmissing <-  which(allwindstrand$chrom == uniquechrom &
                                    allwindstrand$transcript == uniquetrans &
                                    allwindstrand$gene == uniquegene &
                                    allwindstrand$window == idxna)
                            if (!isTRUE(all.equal(length(idxmissing), 1)))
                                stop("Problem in retrieving the missing window, this should not happen. Contact the developper.")
                            
                            ## Below the bedgraph information columns are set to NA. These columns will be removed later
                            ## The score is set to NA since it is a missing value resulting from removing black list and low mappability (keeping high mappability)
                            ## Filling the other columns with the line retrieved in allwindstrand
                            windstrandrow <- allwindstrand[idxmissing, ]
                            resmissing <- data.frame(chrom = windstrandrow$chrom,
                                                start.bg = NA, end.bg = NA, width.bg = NA, strand.bg = "*", ## Set the bedgraph info
                                                score.bg = NA, ## Set the score to NA to keep track of missing values
                                                biotype.window = windstrandrow$biotype, start.window = windstrandrow$start,
                                                end.window = windstrandrow$end, transcript = windstrandrow$transcript, gene = windstrandrow$gene,
                                                strand.window = windstrandrow$strand, window = windstrandrow$window, coord = windstrandrow$coord)
                            return(resmissing)

                        }, allwindstrand)
                        missingrowsdf <- do.call("rbind", missingrowslist)
                        currenttrans <- rbind(currenttrans, missingrowsdf)
                        currenttrans <- currenttrans[order(currenttrans$coord), ]
                    }

                    score.bg <- currenttrans$score.bg
                    currenttrans <- currenttrans[, -grep(".bg", colnames(currenttrans))]
                    currenttrans <- cbind(currenttrans, score.bg)
                    colnames(currenttrans)[which(colnames(currenttrans) == "score.bg")] <- paste0(currentname, ".score")
                    return(currenttrans)
                }, windsize, allwindstrand, currentname)

            }, allwindstrand, currentname, resblack, maptracktib)
        }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, verbose), SIMPLIFY = FALSE,
        mc.cores = nbcpubg)
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                

                ## Setting missing frames to NA
               
                    

                            # misschrom <- currenttrans$chrom[1]
                            # misstrans <- currenttrans$transcript.anno[1]
                            # missgene <- currenttrans$gene.anno[1]

                            # idxmiss <- which(resmap$chrom == misschrom &
                            #     resmap$transcript.anno == misstrans &
                            #     resmap$gene.anno == missgene)
                            
                            # !!!!!!!!! RETRIEVE CURRENTNA - 1 EXCEPT IF CURRENTNA IS 1
                            # !!!!!!!!!!! THEN MODIFY CURRENTTRANS WITH <<-
                            # !!!!!!!!!! CHANGE LOOP TO INVISIBLE
                            # if (!isTRUE(all.equal(length(idxmiss), 1)))
                            #     stop("idxmiss should be unique in retrieveandfilterfrombg. Contact the developper.")
                            
                            # resrow <- resmap[idxmiss, ]
                            # resrow["score"] <- NA

                            
                            if (isTRUE(all.equal(currentna, 1))) {} else {}
                                
                            resrow <- data.frame(chrom = currenttrans$chrom[1]
                            misstrans <- currenttrans$transcript.anno[1]
                            missgene <- currenttrans$gene.anno[1]



               })
                
            }, windsize, resmap, allwindtib)

            !!
            
            return(resmap)
            }, allwindtib)

        return(resmap)

    })

    return(bedgraphlist)
}



.computewmeanvec <- function(dupframenbvec, df, expname, colscore) {
    wmeanvec <- sapply(dupframenbvec, function(namedup, df, expname, colscore) {

        ## Selecting all rows having a duplicated frame found at index idx
        allframedf <- df[which(df$window == namedup), ]
        if (isTRUE(all.equal(nrow(allframedf), 1)))
            stop("There should be more than one frame selected")

        ## Testing that the coord of the window is the same for all scores
        ## selected (this should not give an error)
        if (!isTRUE(all.equal(length(unique(allframedf[, "start"])), 1)) ||
            !isTRUE(all.equal(length(unique(allframedf[, "end"])), 1)))
                stop("The size of the window is not unique for the frame rows ",
                    "selected, this should not happen, contact the developper.")

        ## Retrieving the coordinates and the size of the transcript
        windowstart <- allframedf[1, "start"]
        windowend <- allframedf[1, "end"]

        ## Retrieve the nb of overlapping nt for each score
        overntvec <- apply(allframedf, 1,
            function(x, expname, windowstart, windowend) {
                nt <- seq(from = x[paste0(expname, "start")],
                    to = x[paste0(expname, "end")], by = 1)
                overnt <- length(which(nt >= windowstart & nt <= windowend))
                return(overnt)
            }, expname, windowstart, windowend)

        ## Computing weighted mean
        wmean <- weighted.mean(allframedf[, colscore], overntvec)
        return(wmean)
    }, df, expname, colscore)
    return(wmeanvec)
}





        # tab=bgscorebytrans[[1]]; nametrs=names(idxbgscorebytrans)[1]
        # annogr=allwindowsgr;bggr=currentgr;strd=currentstrand;
        # expname=currentname
        dfwmeanbytranslist <- mcmapply(function(tab, nametrs, annogr, bggr,
            strd, expname, windsize) {

            ## Building the complete data.frame and identifying duplicated
            ## frames
            df <- .buildtransinfotable(annogr, tab, bggr, expname)
            dupidx <- which(duplicated(df$window))
            colscore <- paste0(expname, "score")

            if (!isTRUE(all.equal(length(dupidx), 0))) {
                dupframenbvec <- unique(df$window[dupidx])
                ## For each duplicated frame
                wmeanvec <- .computewmeanvec(dupframenbvec, df, expname,
                    colscore)

                ## Remove duplicated frames and replace scores by wmean and
                ## adding the coord column
                df <- .replaceframeswithwmean(df, dupidx, windsize, nametrs,
                    dupframenbvec, colscore, strd, wmeanvec)
            }

            return(df)
        }, idxbgscorebytrans, names(idxbgscorebytrans),
        MoreArgs = list(allwindowsgr, currentgr, currentstrand, currentname,
        windsize), SIMPLIFY = FALSE, mc.cores = nbcputrans)

        return(dfwmeanbytranslist)
}



        
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



!!!!!!!!!!!!!!!!! CODE FOR KEEPING ON HIGH MAPPABILITY
        resmaplist <- lapply(chromvec, function(currentchrom) {

            if (verbose) message("\t\t\t over ", currentchrom)
            ## Keeping scores on high mappability track
            if (verbose) message("\t Keeping scores on high mappability track")
            resmap <-  tryCatch({
                valr::bed_intersect(resblack,
                maptracktib  %>% dplyr::filter(chrom == currentchrom), # nolint
                suffix = c("", "maphigh"))
            }, error = function(e) {
                message(e)
                return(NA)
            })
            return(resmap)})



!!!!!!!!!!!!!!!
