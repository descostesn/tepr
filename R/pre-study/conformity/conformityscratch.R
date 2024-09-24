library("rtracklayer")
library("GenomicRanges")
library("tibble")
library("dplyr")
library("valr")
library("utils")
library("parallel")
library("tidyr")
library("tidyselect")

##################
# PARAMETERS
##################


#### preprocessing

bgvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score" # nolint
allbgnicpath <- "/g/romebioinfo/tmp/preprocessing/backup/bedgraphwmeanlist.rds"
allwindowspath <- "/g/romebioinfo/tmp/preprocessing/allwindowsbed.rds"
completedfpath <- "/g/romebioinfo/tmp/preprocessing/backup/completeframedf.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv" # nolint
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint

nbcpubg <- 8
nbcputrans <- 20
windsize <- 200

#### downstream

vicbigtsvpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/dTAG_Cugusi_stranded_20230810.tsv" # nolint
expdfpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
nicbigtsvpath <- "/g/romebioinfo/tmp/comparewithscratch/finaltab.rds"
oldnictsvpath <- "/g/romebioinfo/tmp/preprocessing/backup/completeframedf.rds"
expthres <- 0.1
outputfolder <- "/g/romebioinfo/tmp/comparewithscratch"



##################
#FUNCTIONS - preprocessing
##################

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
#FUNCTIONS - downstream
##################


averageandfilterexprs <- function(expdf, alldf, expthres, verbose = FALSE) { # nolint

    ## Adding column names to alldf
    infocolnames <- c("biotype", "chr", "coor1", "coor2", "transcript",
        "gene", "strand", "window", "id")
    expcolnames <- unlist(apply(expdf, 1, function(x) {
        res <- paste0(x["condition"], "_rep", x["replicate"], ".", x["strand"])
        return(c(res, paste(res, "score", sep = "_")))
    }, simplify = FALSE))
    colnames(alldf) <- c(infocolnames, expcolnames)

    scorecolvec <- expcolnames[grep("_score", expcolnames)]

    ## Calculate the average expression per transcript (over each frame)
    if (verbose) message("\t Calculating average expression per transcript") # nolint
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint

    ## Remove a line if it contains only values < expthres (separating strands)
    if (verbose) message("\t Removing lines with values < expthres") # nolint
    dfstrandlist <- mapply(function(strandname, directname, dfbytrans,
        expthres) {
            if ((isTRUE(all.equal(strandname, "plus")) &&
                isTRUE(all.equal(directname, "reverse"))) ||
                (isTRUE(all.equal(strandname, "minus")) &&
                isTRUE(all.equal(directname, "forward"))))
                    stop("Strand and direction do not match, contact the ",
                        "developer")
            dfstrand <- dfbytranscript %>%
                dplyr::filter(strand == strandname) %>% # nolint
                dplyr::select(gene, transcript, strand, # nolint
                tidyselect::contains(directname))  %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ !is.na(.))) %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ . > expthres))
            return(dfstrand)
        }, unique(dfbytranscript$strand), unique(expdf$strand),
            MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)

    exptranstab <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>%
        dplyr::arrange(transcript) %>% dplyr::pull(transcript) # nolint

    return(list(maintable = alldf, exptranlist = exptranstab))
}



.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname, # nolint
                " should be unique, contact the developer.") # nolint
}

.coordandfilter <- function(str, transtable, nbrows) { # nolint

  if (isTRUE(all.equal(str, "minus"))) {
    coordvec <- seq(from = nbrows, to = 1, by = -1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable[order(coordvec), ]
    transtable <- transtable %>% dplyr::select(!dplyr::matches("plus"))
  } else {
    coordvec <- seq(from = 1, to = nbrows, by = 1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable %>% dplyr::select(!dplyr::matches("minus"))
  }
  return(transtable)
}

.extractstr <- function(transtable) {

    str <- as.character(unique(transtable$strand))
    .checkunique(str, "str")
    if (isTRUE(all.equal(str, "+"))) {
        str <- "plus"
    } else if (isTRUE(all.equal(str, "-"))) {
        str <- "minus"
    } else {
        stop("In .computeecdf, strand is neither plus or minus. This ",
            "should not happen. Contact the developer.")
    }
    return(str)
}

.computeecdf <- function(transtable, expdf, rounding, nbrows) { # nolint

        ## Retrieving keyword plus or minus
        str <- .extractstr(transtable)

        ## Create the coordinate column and select scores having the righ
        ## orientation
        transtable <- .coordandfilter(str, transtable, nbrows)
        colnamevec <- colnames(transtable)
        colscorevec <- colnamevec[grep("_score", colnamevec)]

        ## Filling the NA of the score columns of the right strand with
        ## tidyr::fill in the downup direction
        transtable <- transtable %>% tidyr::fill(tidyr::contains("score"),
           .direction = "downup")

        ## Computing ecdf
        suppressWarnings(dflong <- transtable %>%
            tidyr::gather(key = "variable", value = "value", colscorevec))
        dflong[, "value_round"] <- round(dflong$value * rounding)
        ecdflist <- lapply(unique(dflong$variable), function(currentvar) {
            dfsubset <- subset(dflong, subset = variable == currentvar) # nolint
            dfexpanded <- dfsubset[rep(seq_len(nrow(dfsubset)),
                dfsubset$value_round), ]
            funecdf <- ecdf(dfexpanded[, "coord"])
            dfsubset$Fx <- funecdf(dfsubset$coord)
            return(dfsubset)
        })
        resecdf <- dplyr::bind_rows(ecdflist)

        ## Shrink the results back to the transtable keeping ecdf columns
        res <- resecdf %>% tidyr::pivot_wider(.,
            names_from = "variable",
            values_from = c("value", "value_round", "Fx")) %>%
            dplyr::select(., -tidyselect::contains("value_round"))

        ## Removing strand from column names
        res <- res %>% dplyr::rename_with(~gsub(paste0(".", str), "", .),
                tidyselect::contains(paste0(".", str)))

        return(res)
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1, # nolint
  verbose = FALSE) {

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]
    maincolnamevec <- colnames(maintable)

    ## Filtering the main table to keep only the expressed transcripts
    if (verbose) message("\t Filtering to keep only the expressed transcripts") # nolint
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    if (isTRUE(all.equal(length(idxnoexpr), 0)))
      warning("All the transcripts are expressed", immediate. = TRUE) # nolint
    else
      maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    if (verbose) message("\t Splitting the table by each transcript") # nolint
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        rounding, nbrows, maincolnamevec) {

        res <- .computeecdf(transtable, expdf, rounding, nbrows)
        return(res)
    }, expdf, rounding, nbrows, maincolnamevec, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(list(concatdf, nbrows))
}


.condcolidx <- function(currentcond, df) {
    idxcond <- grep(currentcond, colnames(df))
    if (isTRUE(all.equal(length(idxcond), 0)))
        stop("Problem in function createmeandiff, condition not found in ",
                "column names. Contact the developer.")
    return(idxcond)
}

.idxscorefx <- function(df, idxcond) {
    idxcondfx <- grep("Fx", colnames(df[idxcond]))
    idxcondval <- grep("value_", colnames(df[idxcond]))
    if (isTRUE(all.equal(length(idxcondfx), 0)) ||
        isTRUE(all.equal(length(idxcondval), 0)))
        stop("Problem in function createmeandiff, column Fx or val not found ",
            "in column names. Contact the developer.")
    idxcondlist <- list(value = idxcond[idxcondval],
            Fx = idxcond[idxcondfx])
    return(idxcondlist)
}

.creatematdiff <- function(condvec, resmean) {

  categoryvec <- c("value", "Fx")
  matdifflist <- lapply(categoryvec, function(currentcat, condvec, resmean) {
    meancolnames <- paste("mean", currentcat, condvec, sep = "_")

    ## Generating all combinations of elements (combn not good)
    idxvec <- seq_len(length(condvec))
    matidx <- matrix(c(idxvec, rev(idxvec)), ncol = 2)

    ## Generating differences of columns
    difflist <- apply(matidx, 2, function(idxvec, meancolnames, resmean,
        currentcat, condvec) {
          ## The function rowDiffs of the package matrixStats substracts the second argument to the first one. To respect the code just above, # nolint
          ## The indexes must be inverted with rev: meancolnames[rev(idxvec)]] -> for instance, given the two columns "mean_value_HS" and "mean_value_ctrl" # nolint
          ## as input, the function rowDiffs will do the subtraction "mean_value_ctrl" - "mean_value_HS" # nolint
          res <- matrixStats::rowDiffs(as.matrix(
            resmean[, meancolnames[rev(idxvec)]]))
          colnamestr <- paste("Diff", paste0("mean", currentcat),
            paste(condvec[idxvec], collapse = "_"), sep = "_")
          res <- as.vector(res)
          attr(res, "name") <- colnamestr
          return(res)
      }, meancolnames, resmean, currentcat, condvec, simplify = FALSE)

      ## Combining vectors into a matrix and defining col names
      diffmat <- do.call("cbind", difflist)
      colnames(diffmat) <- sapply(difflist, function(x) attributes(x)$name)
      return(diffmat)
  }, condvec, resmean)

  ## Building a matrix from the diff on values and Fx
  matdiff <- do.call("cbind", matdifflist)
  return(matdiff)
}

.meandiffscorefx <- function(idxcondlist, df, nbwindows, currentcond,
    colnamevec, verbose) {

        meandifflist <- mapply(function(idxvalvec, idxname, df, nbwindows,
            currentcond, colnamevec, verbose) {
            if (verbose) {
              message("\t Calculating average and difference between ",
                "replicates for columns '", idxname, "' of ", currentcond)
              if (isTRUE(all.equal(length(idxvalvec), 1)))
                warning("Only one replicate, copy scores to mean columns",
                  immediate. = TRUE)
            }

            ## Calculating the column of mean scores for currentcond
            ## The result is a data.frame made of a single column
            if (length(idxvalvec) >= 2) {
                meandf <- data.frame(rowMeans(df[, idxvalvec], na.rm = FALSE))
            } else {
                meandf <- as.data.frame(df[, idxvalvec])
            }
            colnames(meandf) <- paste0("mean_", idxname, "_", currentcond)

            if (isTRUE(all.equal(idxname, "Fx"))) {
                diffres <- meandf - df$coord / nbwindows
                colnames(diffres) <- paste0("diff_", idxname, "_", currentcond)
                res <- cbind(meandf, diffres)
            } else {
                res <- meandf
            }
            return(res)
        }, idxcondlist, names(idxcondlist), MoreArgs = list(df, nbwindows,
            currentcond, colnamevec, verbose), SIMPLIFY = FALSE)

        return(meandifflist)
}

createmeandiff <- function(resultsecdf, expdf, nbwindows, verbose = FALSE) {

    ## for each condition, creates three columns:
    ##   - "mean_value_ctrl", "mean_Fx_ctrl", "diff_Fx_ctrl"
    ##   - "mean_value_HS", "mean_Fx_HS", "diff_Fx_HS"
    condvec <- unique(expdf$condition)
    rescondlist <- lapply(condvec, function(currentcond, df, nbwindows,
      verbose) {

        if (verbose) message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- .condcolidx(currentcond, df)
        ## Separating idx of column names by scores and Fx
        idxcondlist <- .idxscorefx(df, idxcond)

        ## The difference is used to calculate the AUC later on
        colnamevec <- colnames(df)
        meandifflist <- .meandiffscorefx(idxcondlist, df, nbwindows,
          currentcond, colnamevec, verbose)
        names(meandifflist) <- NULL

        meandiffres <- do.call("cbind", meandifflist)
        return(meandiffres)
    }, resultsecdf, nbwindows, verbose)
    resmean <- do.call("cbind", rescondlist)

#!!!!!!!!!!!!!  CODE TO REMOVE FOR THE FINAL FUNCTION
        # res <- cbind(resultsecdf, resmean)
        # return(res)

    ## Computing all differences on mean columns
    if (verbose) message("Commputing all differences on mean columns")
    matdiff <- .creatematdiff(condvec, resmean)

    res <- cbind(resmean, matdiff)
    if (!isTRUE(all.equal(nrow(resultsecdf), nrow(res))))
        stop("The results of mean and diff should have the same number of ",
            "rows than resultsecdf, contact the developer")

    return(cbind(resultsecdf, res))
}


.returninfodf <- function(transtab, nbwindows) { # nolint

    infodf <- transtab  %>%
        dplyr::filter(window == round(nbwindows / 2))  %>%
        dplyr::mutate(window_size = abs(coor2 - coor1), .keep = "all") %>% # nolint
        dplyr::select("transcript", "gene", "strand", "window_size") %>%
        dplyr::distinct()

    return(infodf)
}

.checkempty <- function(idx, namestr) {

    if (isTRUE(all.equal(length(idx), 0)))
        stop("Your condition ", namestr, " was not found in the ",
            "experiment table expdf. Please verify")
}

.dauc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1,
    controlcondname = "ctrl", stresscondname = "HS", dontcompare = NULL) {

    condvec <- unique(expdf$condition)
    resdflist <- mclapply(bytranslist, function(transtab, condvec,
        controlcondname, stresscondname) {

        ## Retrieve the column names for each comparison
        idxctrl <- grep(controlcondname, condvec)
        .checkempty(idxctrl, controlcondname)
        idxstress <- grep(stresscondname, condvec)
        .checkempty(idxstress, stresscondname)
        name1 <- paste0("mean_Fx_", condvec[idxctrl])  # nolint
        name2 <- paste0("mean_Fx_", condvec[idxstress])  # nolint
        diffname <- paste0("Diff_meanFx_", condvec[idxstress], "_",  # nolint
          condvec[idxctrl])

        ## Perform a kolmogorov-smirnoff test between the two columns
        resks <- suppressWarnings(ks.test(transtab[, name1], transtab[, name2]))

        ## Calculate the area under the curve of the difference of means
        ## -> delta AUC
        deltadauc <- pracma::trapz(transtab[, "coord"], transtab[, diffname])
        ## Retrieve the p-value
        pvaldeltadaucks <- resks$p.value
        ## The KS test statistic is defined as the maximum value of the
        ## difference between A and B’s cumulative distribution functions (CDF)
        statdeltadaucks <- resks$statistic

        ## Build a one line data.frame with the proper col names
        ksdeltadaucdf <- data.frame(deltadauc, pvaldeltadaucks, statdeltadaucks)
        prefixvec <- c("dAUC", "p_dAUC", "D_dAUC")
        colnames(ksdeltadaucdf) <- paste(prefixvec, diffname, sep = "_")

        ## Retrieving transcript information
        infodf <- .returninfodf(transtab, nbwindows)

        ## Combining the two df as result
        resdf <- cbind(infodf, ksdeltadaucdf)
        return(resdf)
    }, condvec, controlcondname, stresscondname, mc.cores = nbcpu)

    resdf <- do.call("rbind", resdflist)

    ## Correct p-values using FDR
    idx <- grep("p_dAUC", colnames(resdf))
    fdrvec <- p.adjust(resdf[, idx], method = "fdr")

    resdf <- cbind(resdf, fdrvec)
    colnamevec <- colnames(resdf)
    idxfdr <- which(colnamevec == "fdrvec")
    colnames(resdf)[idxfdr] <- paste0("adjFDR_", colnamevec[idx])  # nolint
    return(resdf)
}

.buildaucdf <- function(transtab, difffxname, resks, meanvalname,
  currentcond, nbwindows) {
    auc <- pracma::trapz(transtab[, "coord"], transtab[, difffxname])
    pvalaucks <- resks$p.value
    stataucks <- resks$statistic
    meanvaluefull <- mean(transtab[, meanvalname])
    aucdf <- data.frame(auc, pvalaucks, stataucks, meanvaluefull)
    prefixvec <- c("AUC", "p_AUC", "D_AUC", "MeanValueFull")
    colnames(aucdf) <- paste(prefixvec, currentcond, sep = "_")
    rownames(aucdf) <- NULL
    transinfo <- .returninfodf(transtab, nbwindows)
    res <- cbind(transinfo, aucdf)
    return(res)
}

.auc_allconditions <- function(bytranslist, expdf, nbwindows, nbcpu = 1) {

  cumulative <- seq(1, nbwindows) / nbwindows
  condvec <- unique(expdf$condition)

  resdflist <- mclapply(bytranslist, function(transtab, condvec, cumulative,
    nbwindows) {
      ## Computing AUC, pval, and stat for each condition
      resauclist <- lapply(condvec, function(currentcond, transtab,
        cumulative, nbwindows) {
          ## Definition of column names
          difffxname <- paste0("diff_Fx_", currentcond) # nolint
          meanvalname <- paste0("mean_value_", currentcond) # nolint
          meanfxname <- paste0("mean_Fx_", currentcond) # nolint

          ## Perform a kolmogorov-smirnoff test between mean_Fx and cum.density
          resks <- suppressWarnings(ks.test(transtab[, meanfxname], cumulative))
          ## Build data.frame with auc information for the current transcript
          aucdf <- .buildaucdf(transtab, difffxname, resks, meanvalname,
            currentcond, nbwindows)
          return(aucdf)
      }, transtab, cumulative, nbwindows)
      aucdf <- do.call("cbind", resauclist)
      return(aucdf)
  }, condvec, cumulative, nbwindows, mc.cores = nbcpu)

  aucallconditions <- do.call("rbind", resdflist)
  idxdup <- which(duplicated(colnames(aucallconditions)))
  aucallconditions <- aucallconditions[, -idxdup]

  ## Correcting p-val with FDR
  idxpvalvec <- grep("p_AUC", colnames(aucallconditions))
  fdrlist <- lapply(idxpvalvec, function(idxpval, tab) {
    return(p.adjust(tab[, idxpval], method = "fdr"))
  }, aucallconditions)
  fdrdf <- do.call("cbind", fdrlist)
  colnames(fdrdf) <- paste0("adjFDR_", colnames(aucallconditions)[idxpvalvec]) # nolint

  aucallconditions <- cbind(aucallconditions, fdrdf)
  return(aucallconditions)
}

allauc <- function(bytranslistmean, expdf, nbwindows, nbcputrans,
  dontcompare = NULL, controlcondname = "ctrl", stresscondname = "HS",
  verbose = TRUE) {

    if (isTRUE(all.equal(length(unique(expdf$condition)), 2))) {
        if (verbose) message("\t Computing the differences (d or delta) of AUC")
        start_time <- Sys.time()
        daucallcond <- .dauc_allconditions(bytranslistmean, expdf, nbwindows,
          nbcputrans, controlcondname, stresscondname)
        end_time <- Sys.time()
        if (verbose) message("\t\t ## Analysis performed in: ",
          end_time - start_time) # nolint
    } else {
        warning("dAUC not performed, only one condition submitted.")
    }

## !!!!!!!!!!!!!!! TO REMOVE FROM FINAL CODE
##        return(daucallcond)


    ## Calculate the Area Under Curve (AUC), All conditions vs y=x
    ## Calculate Mean Value over the full gene body in All conditions.
    if (verbose) message("\t Computing the Area Under Curve (AUC)")
    start_time <- Sys.time()
    aucallcond <- .auc_allconditions(bytranslistmean, expdf, nbwindows,
      nbcpu = nbcputrans)
    end_time <- Sys.time()
    if (verbose) message("\t\t ## Analysis performed in: ",
      end_time - start_time) # nolint
    #saveRDS(aucallcond, "/g/romebioinfo/tmp/downstream/aucallcond.rds") # nolint

## !!!!!!!!!!!!!!! TO REMOVE FROM FINAL CODE
##        return(aucallcond)

    ## Merging the two tables by transcript
    if (verbose) message("Merging results")
    allauc <- merge(aucallcond, daucallcond,
      by = c("gene", "transcript", "strand", "window_size"))
    return(allauc)
}




##################
# MAIN
##################

#########################################
# PART 1: PREPROCESSING
#########################################


## Reading all windows bed
allwindowsbed <- readRDS(allwindowspath)

## Reading exptab, black list, and maptrack
exptab <- read.csv(exptabpath, header = TRUE)
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
blacklistbed <- read.delim(blacklistshpath, header = FALSE)
maptrackbed <- read.delim(maptrackpath, header = FALSE)

bedgraphlistwmean <- retrieveandfilterfrombg(exptab, blacklistbed, maptrackbed,
    nbcputrans, allwindowsbed, expnamevec, windsize)

message("Merging results of each bedgraph into a single table")
finaltab <- createtablescores(bedgraphlistwmean, nbcpubg)
saveRDS(finaltab, file = file.path(outputfolder, "finaltab.rds"))



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


#########################################
# PART 2: DOWNSTREAM
#########################################


## This code tests the functions of downstream.R with the input table of
## victor: /g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/dTAG_Cugusi_stranded_20230810.tsv # nolint
bigtsv <- read.table(vicbigtsvpath, header = FALSE)
expdf <- read.csv(expdfpath, header = TRUE)


####
#### averageandfilterexprs
####

niccode_allexprsdfsvic <- averageandfilterexprs(expdf, bigtsv, expthres,
    verbose = TRUE)
viccode_allexprsdfsvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/results_main_table.rds") # nolint

if (isTRUE(all.equal(niccode_allexprsdfsvic[[1]], viccode_allexprsdfsvic[[1]])))
    message("table is equal after averageandfilterexprs")

if (isTRUE(all.equal(niccode_allexprsdfsvic[[2]], viccode_allexprsdfsvic[[2]])))
    message("transcript list is equal after averageandfilterexprs")


####
#### genesECDF
####

niccode_resecdfvic <- genesECDF(niccode_allexprsdfsvic, expdf, nbcpu = nbcputrans, verbose = TRUE) # nolint
nbwindows <- niccode_resecdfvic[[2]]
niccode_resecdfvic <- niccode_resecdfvic[[1]]

## Reading the result of ecdf that contains the column coord that is present in
## the input table of nic
viccode_resecdfvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/cugusi2023_ECDFScores_10_200.tsv" # nolint
viccode_resecdfvic <- read.table(viccode_resecdfvicpath, sep = "\t", header = TRUE) # nolint

if (isTRUE(all.equal(as.data.frame(niccode_resecdfvic), viccode_resecdfvic)))
    message("genesECDF is consistent")


####
#### createmeandiff
####

## IMPORTANT: For the sake of comparison with the code of vic, only the first
## part of createmeandiff was executed by adding the following lines after
## resmean:
##        res <- cbind(resultsecdf, resmean)
##        return(res)

viccode_dfmeanvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/concat_dfFX_res.rds") # nolint
niccode_dfmeanvic <- createmeandiff(niccode_resecdfvic, expdf, nbwindows)

if (isTRUE(all.equal(viccode_dfmeanvic, niccode_dfmeanvic)))
    message("consistancy after dfmean")

## IMPORTANT: Now the whole function is executed (above lines are commented) to
## compute the differences of means
viccode_dfmeandiffvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/concat_Diff_mean_res.rds") # nolint
niccode_dfmeandiffvic <- createmeandiff(niccode_resecdfvic, expdf, nbwindows)

## Change the 'V' of 'value' to lower case in vic table
colnames(viccode_dfmeandiffvic)[which(colnames(viccode_dfmeandiffvic) == "Diff_meanValue_ctrl_HS")] <- "Diff_meanvalue_ctrl_HS" # nolint
colnames(viccode_dfmeandiffvic)[which(colnames(viccode_dfmeandiffvic) == "Diff_meanValue_HS_ctrl")] <- "Diff_meanvalue_HS_ctrl" # nolint

if (isTRUE(all.equal(viccode_dfmeandiffvic, niccode_dfmeandiffvic)))
    message("consistancy after mean differences")


####
#### dAUC
####


bytranslistmean <- split(niccode_dfmeandiffvic,
    factor(niccode_dfmeandiffvic$transcript))

## IMPORTANT: Only the delta auc (first part) of the function allauc is executed
## by adding return(daucallcond) in the function. This should be removed later.
niccode_daucdfvic <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)
viccode_daucdfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/dAUC_allcondi_res.rds") # nolint

# For comparison only
rownames(niccode_daucdfvic) <- NULL
names(viccode_daucdfvic$D_dAUC_Diff_meanFx_HS_ctrl) <- NULL

if (isTRUE(all.equal(viccode_daucdfvic, niccode_daucdfvic)))
    message("consistancy after dAUC")

####
#### AUC
####

## IMPORTANT: Only the second part of allauc is executed here by skipping the
## code when pasting in the terminal and adding return(aucallcond)
niccode_aucdfvic <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)
viccode_aucdfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/AUC_allcondi_res.rds") # nolint

# For comparison only
rownames(niccode_aucdfvic) <- rownames(viccode_aucdfvic) <- NULL
names(viccode_aucdfvic$D_AUC_ctrl) <- names(viccode_aucdfvic$D_AUC_HS) <- NULL

if (isTRUE(all.equal(viccode_aucdfvic, niccode_aucdfvic)))
    message("consistancy after AUC")


####
#### countNA
####

!!!!!!!!!!!!!!!!!!!!!!!!!!!
countna <- function(allexprsdfs, expdf, nbcpu, verbose = FALSE) {

  maintable <- allexprsdfs[[1]]
  scorecolvec <- grep("_score", colnames(maintable), value = TRUE)
  condvec <- unique(expdf$condition)

  ## Splitting the table by each transcript
  if (verbose) message("\t Splitting the table by each transcript") # nolint
  transdflist <- split(maintable, factor(maintable$transcript))

  ## For each transcript
  nabytranslist <- parallel::mclapply(transdflist,
    function(transtable, scorecolvec, condvec) {
        ## Filters the score columns according to the strand of the transcript
        str <- .extractstr(transtable)
        colnamestr <- scorecolvec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]

        ## Counting NA for each condition (c=condition, m=matrix, n=colnames)
        res <- sapply(condvec, function(c, m, n) {
          idxcol <- grep(c, n)
          if (!isTRUE(all.equal(length(idxcol), 1))) {
            return(length(which(apply(m[, idxcol], 1,
              function(x) all(is.na(x))))))
          } else { ## Only one replicate
            return(length(which(is.na(m[, idxcol]))))
          }
        }, scoremat, colnamestr)

        ## Retrieving total NA and transcript info
        if (!isTRUE(all.equal(res[1], res[2])))
            stop("Number of NA is different between conditions. This should ",
                "not happen. Contact the developer.")
        ## I drop the other NA columns because it is the same value for all the
        ## conditions (NA depends on blacklist and unmmapable region)
        countna <- res[1]
        info <- data.frame(gene = unique(transtable$gene),
          transcript = unique(transtable$transcript),
          strand = unique(transtable$strand))
        return(cbind(info, Count_NA = countna))
    }, scorecolvec, condvec, mc.cores = nbcpu)

  return(do.call("rbind", nabytranslist))
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
niccode_countnavic <- countna(niccode_allexprsdfsvic, expdf, nbcputrans)

viccode_countnavic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/count_NA_res.rds") # nolint