.wmeanvec <- function(dupframenbvec, currenttrans) {

    wmeanvec <- sapply(dupframenbvec, function(nbdup, currenttrans) {

        ## Selecting all rows having a window equal to nbdup
        allframedf <- currenttrans[which(currenttrans$window.window == nbdup), ]

        ## Retrieving start and end of the window to calculate nb of nt
        windowstart <- unique(allframedf$start.window)
        windowend <- unique(allframedf$end.window)
        if (!isTRUE(all.equal(length(windowstart), 1)) ||
            !isTRUE(all.equal(length(windowend), 1)))
            stop("The size of the window is not unique for the frame rows ",
                "selected, this should not happen, contact the developper.")

        ## Retrieve the nb of overlapping nt for each score
        overntvec <- apply(allframedf, 1, function(x, windowstart, windowend) {
            nt <- seq(from = x["start"], to = x["end"], by = 1)
            overnt <- length(which(nt >= windowstart & nt <= windowend))
            return(overnt)
        }, windowstart, windowend)

        ## Computing weighted mean
        allscores <- as.data.frame(allframedf[, "score"])[[1]]
        wmean <- weighted.mean(allscores, overntvec)
        return(wmean)
    }, currenttrans)

    return(wmeanvec)
}

.removeblackandlowmap <- function(currenttrans, blacklisttib, idxscore,
    maptracktib) {

        ## Set scores overlapping black list to NA
        resblack <- valr::bed_intersect(currenttrans, blacklisttib)
        if (!isTRUE(all.equal(nrow(resblack), 0))) {
            strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                currenttrans$end, sep = "-")
            strblack <- paste(resblack$chrom, resblack$start.x, resblack$end.x,
                sep = "-")
            idxblack <- as.vector(na.omit(unique(match(strtransvec, strblack))))
            if (isTRUE(all.equal(length(idxblack), 0)))
                stop("Problem in setting scores overlapping black list to NA.",
                    " This should not happen. Contact the developer.")
            currenttrans[idxblack, idxscore] <- NA
        }

        ## Set scores NOT overlapping high map to NA (i.e. scores overlapping
        ## low mappability intervals)
        resmap <- valr::bed_intersect(currenttrans, maptracktib)
        if (!isTRUE(all.equal(nrow(resmap), 0))) {
            ## Compute strtransvec only if no overlap with black list was found
            if (isTRUE(all.equal(nrow(resblack), 0)))
                strtransvec <- paste(currenttrans$chrom, currenttrans$start,
                    currenttrans$end, sep = "-")
            strmap <- paste(resmap$chrom, resmap$start.x, resmap$end.x,
                sep = "-")
            idxmap <- match(strtransvec, strmap)
            idxtosetna <- which(is.na(idxmap)) ## NOT overlapping high map
            currenttrans[idxtosetna, idxscore] <- NA
        }
        return(currenttrans)
}

.meanblackhighbytrans <- function(bgscorebytrans, windsize, currentname,
    currentchrom, blacklisttib, maptracktib, saveobjectpath, nbcputrans,
    reload, verbose) {

        currentobj <- file.path(saveobjectpath, paste0(currentname,
            "-", currentchrom, "-translist.rds"))
        if (!reload || !file.exists(currentobj)) {
            bytranslist <- parallel::mclapply(bgscorebytrans,
                function(currenttrans, windsize, currentname, blacklisttib,
                    maptracktib) {

                        ## Reordering rows according to window number
                        idxwind <- order(currenttrans$window.window)
                        currenttrans <- currenttrans[idxwind, ]

                        ## Identifying duplicated windows that will be used to
                        ## compute a weighted mean.
                        currenttrans <- .dupidx(currenttrans, windsize)

                        ## Remove columns corresponding to bedgraph, move score
                        ## column at the end, remove the .window suffix from
                        ## column names, Add exp name prefix to column score
                        res <- .formatcurrenttranscols(currenttrans,
                            currentname)
                        currenttrans <- res[[1]]
                        idxscore <- res[[2]]

                        ## Set scores overlapping black list and low map to NA
                        idxchrom <- which(maptracktib$chrom == unique(
                            currenttrans$chrom))
                        maptracktibchrom <- maptracktib[idxchrom, ]
                        currenttrans <- .removeblackandlowmap(currenttrans,
                            blacklisttib, idxscore, maptracktibchrom)

                        return(currenttrans)

            }, windsize, currentname, blacklisttib, maptracktib,
                mc.cores = nbcputrans)
            if (!is.na(saveobjectpath)) {
                if (verbose) message("\t\t\t Saving ", currentobj)
                saveRDS(bytranslist, file = currentobj)
            }
        } else {
            if (verbose) message("\t\t\t Loading ", currentobj)
            bytranslist <- readRDS(currentobj)
        }

        return(bytranslist)
}

.retrieveandfilterfrombg <- function(exptab, blacklisttib, maptracktib, # nolint
    nbcputrans, allwindchromtib, expnamevec, windsize, currentchrom,
    chromlength, saveobjectpath, showtime, reload, tmpfold, verbose) {

        ## Looping on each experiment bg file
        if (verbose) message("\t\t For each bedgraph file") # nolint
        invisible(mapply(function(currentpath, currentname,
            currentstrand, currentcond, currentrep, currentdirection,
            allwindchromtib, blacklisttib, maptracktib, windsize, currentchrom,
            chromlength, nbcputrans, saveobjectpath, verbose, showtime,
            reload, tmpfold) {

            filename <- file.path(tmpfold, paste0(currentname, "-",
                currentchrom, ".tsv"))

            ## Deleting res which is created at the end of the loop before
            ## creating the new one
            if (exists("res")) {
                rm(res)
                invisible(gc())
            }

            if (!reload || !file.exists(filename)) {
                ## Retrieving bedgraph values
                if (verbose) message("\n\t\t Retrieving begraph values for ", # nolint
                    currentname, " on ", currentchrom) # nolint
                valtib <- .retrievebgval(currentpath, currentchrom, chromlength,
                    verbose)

                ## Retrieving scores on annotations of strand
                annoscores <- .retrieveannoscores(currentstrand,
                    allwindchromtib, valtib, verbose)

                ## Splitting the scores by transcript
                if (verbose) message("\t\t Splitting the scores by transcript")
                trsfact <- factor(annoscores$transcript.window)
                bgscorebytrans <- split(annoscores, trsfact)
                rm(trsfact, annoscores)
                invisible(gc())

                ## For each transcript compute the weighted means for each
                ## window. The weight is calculated if a window contains more
                ## than one score
                if (verbose) message("\t\t For each transcript compute the ",
                    "weighted means and set scores overlapping black list and ",
                    "low mappability to NA. It takes a while.")
                if (showtime) start_time_bytranslist <- Sys.time()

                bytranslist <- .meanblackhighbytrans(bgscorebytrans, windsize,
                    currentname, currentchrom, blacklisttib, maptracktib,
                    saveobjectpath, nbcputrans, reload, verbose)

                if (showtime) {
                    end_time_bytranslist <- Sys.time()
                    timing <- end_time_bytranslist - start_time_bytranslist
                    message("\t\t\t ## Features excluded in: ", timing) # nolint
                }

                if (!isTRUE(all.equal(unique(sapply(bytranslist, nrow)),
                    windsize)))
                    stop("All elements of the list should contain ", windsize,
                        " rows. This should not happen. Contact the developer.")
                rm(bgscorebytrans)
                invisible(gc())

                ## Formatting columns and adding rowid column
                res <- .rowidandcols(bytranslist, currentcond, currentrep,
                    currentdirection, verbose)

                ## Saving table to temporary folder
                if (verbose) message("\t\t Saving table to ", filename)
                res <- res[order(res$rowid, res$window), ]
                write.table(res, file = filename, sep = "\t", quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
            } else {
                if (verbose) message("\t\t The file ", filename,
                    " was already computed. Skipping.")
            }

        }, exptab$path, expnamevec, exptab$strand, exptab$condition,
            exptab$replicate, exptab$direction, MoreArgs = list(allwindchromtib,
            blacklisttib, maptracktib, windsize, currentchrom, chromlength,
            nbcputrans, saveobjectpath, verbose, showtime, reload, tmpfold),
            SIMPLIFY = FALSE))
}

.loadbgprocessing <- function(exptab, blacklisttib, maptrackpath, allwindtib,
        windsize, chromtab, nbcputrans, showtime, saveobjectpath, reload,
        tmpfold, verbose) {

            if (verbose) message("Removing scores within black list intervals,",
            " keeping those on high mappability regions, and computing ",
            "weighted means.")
            expnamevec <- paste0(exptab$condition, exptab$replicate,
                exptab$direction)

            ## Loading process on a specific chromosom
            invisible(lapply(GenomeInfoDb::seqnames(chromtab),
                function(currentchrom, chromtab, maptrackpath, showtime,
                saveobjectpath, reload, verbose, exptab, blacklisttib,
                nbcputrans, allwindtib, expnamevec, windsize, tmpfold) {

                    if (showtime) start_bglistwmean <- Sys.time()

                    if (verbose) message("\n\t #--- Processing ", currentchrom)
                    ## Reading the maptrack on a specific chromosomes
                    chromlength <- .retrievechromlength(chromtab, currentchrom)
                    maptracktib <- .retrievemaptrack(maptrackpath, showtime,
                        currentchrom, chromlength, saveobjectpath, reload,
                        verbose)

                    ## Filtering allwindtib on the current chromosome
                    if (verbose) message("\t\t Selecting windows on ",
                        currentchrom)
                    idxchrom <- which(allwindtib$chrom == currentchrom)

                    if (!isTRUE(all.equal(length(idxchrom), 0))) {
                        allwindchromtib <- allwindtib[idxchrom, ]

                        .retrieveandfilterfrombg(exptab, blacklisttib,
                            maptracktib, nbcputrans, allwindchromtib,
                            expnamevec, windsize, currentchrom, chromlength,
                            saveobjectpath, showtime, reload, tmpfold, verbose)

                        if (showtime) {
                            end_bglistwmean <- Sys.time()
                            timing <- end_bglistwmean - start_bglistwmean
                            message("\t\t ## Built ", currentchrom," in: ",
                                timing)
                        }
                    } else {
                        if (verbose) message("\t\t No transcript annotations ",
                            " found on ", currentchrom, ". Skipping.")
                    }

                }, chromtab, maptrackpath, showtime, saveobjectpath, reload,
                    verbose, exptab, blacklisttib, nbcputrans, allwindtib,
                    expnamevec, windsize, tmpfold))
}


# !! DOC TO DO

blacklisthighmap <- function(maptrackpath, blacklistshpath, exptabpath,
    nbcputrans, allwindowsbed, windsize, genomename, saveobjectpath = NA,
    tmpfold = "./tmp", reload = FALSE, showtime = FALSE, verbose = TRUE) {

        if (showtime) start_time_fun <- Sys.time()

        if (!file.exists(tmpfold))
            dir.create(tmpfold, recursive = TRUE)

        ## Retrieving chromosome lengths
        chromtab <- .retrievechrom(genomename, verbose)

        ## Reading the information about experiments
        if (verbose) message("Reading the information about experiments")
        exptab <- read.csv(exptabpath, header = TRUE)

        if (verbose) message("Reading the black list")
        blacklistbed <- read.delim(blacklistshpath, header = FALSE)
        colnames(blacklistbed) <- c("chrom", "start", "end", "type")
        blacklisttib <- tibble::as_tibble(blacklistbed)

        ## Converting allwindowsbed to tib
        colnames(allwindowsbed) <- c("biotype", "chrom", "start", "end",
            "transcript", "gene", "strand", "window")
        allwindtib <- tibble::as_tibble(allwindowsbed)

        ## Cleaning
        rm(blacklistbed, allwindowsbed)
        invisible(gc())

        ## Removing scores within black list intervals, keeping those on high
        ## mappability regions, and computing weighted means.
        .loadbgprocessing(exptab, blacklisttib, maptrackpath, allwindtib,
            windsize, chromtab, nbcputrans, showtime, saveobjectpath, reload,
            tmpfold, verbose)

        if (showtime) {
            end_time_fun <- Sys.time()
            timing <- end_time_fun - start_time_fun
            message("\t\t ## All bedgraphs processed in: ", timing) # nolint
        }
}
