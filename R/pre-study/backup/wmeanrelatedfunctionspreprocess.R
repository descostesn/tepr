## This is a backup of functions that were in preprocessing.R
## new functions have been dev in conformityscratch.R

retrieveandfilterfrombg <- function(exptab, blacklistgr, maptrackgr, nbcpu,
    expnamevec, verbose = TRUE) {

    ## Looping on each experiment bw file
    # currentpath <- exptab$path[1]
    # currentname <- expnamevec[1]
    bedgraphgrlist <- mcmapply(function(currentpath, currentname, blacklistgr,
        maptrackgr, nbcpu, verbose) {

        if (verbose) message("\t Retrieving values for ", currentname)
        valgr <- rtracklayer::import.bedGraph(currentpath)

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

    }, exptab$path, expnamevec, MoreArgs = list(blacklistgr,
        maptrackgr, nbcpu, verbose), mc.cores = nbcpu, SIMPLIFY = FALSE)

    return(bedgraphgrlist)
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

.buildtransinfotable <- function(annogr, tab, bggr, expname) {

    ## Retrieving information about tables
    transcriptvec <- names(annogr)[tab$subjectHits]
    names(annogr) <- NULL
    annodf <- as.data.frame(annogr[tab$subjectHits])
    annodf <- cbind(annodf, transcript = transcriptvec)
    bgdf <- as.data.frame(bggr[tab$queryHits])
    colnames(bgdf) <- paste0(expname, colnames(bgdf))

    ## Building the complete data.frame
    df <- cbind(annodf, bgdf)
    return(df)
}

.replaceframeswithwmean <- function(df, dupidx, windsize, nametrs,
    dupframenbvec, colscore, strd, wmeanvec) {

        ## Remove duplicated frames and replace scores by wmean
        df <- df[-dupidx, ]
        if (!isTRUE(all.equal(nrow(df), windsize)))
            stop("The number of frames should be equal to windsize: ",
                windsize, " for transcript ", nametrs)
        idxscorereplace <- match(dupframenbvec, df$window)
        if (!isTRUE(all.equal(dupframenbvec, df$window[idxscorereplace])))
            stop("Problem in replacing scores by wmean, contact the developer.")
        df[idxscorereplace, colscore] <- wmeanvec

        return(df)
}

summarizebywmean <- function(idxbgscorebytrans, allwindowsgr, currentgr,
    currentstrand, currentname, windsize, nbcputrans) {

        # tab=idxbgscorebytrans[[1]]; nametrs=names(idxbgscorebytrans)[1]
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
