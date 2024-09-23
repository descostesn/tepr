library("rtracklayer")
library("GenomicRanges")
library("tibble")
library("dplyr")
library("valr")
library("utils")
library("parallel")



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

!!!!!!!!!!!!!!

.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname, # nolint
                " should be unique, contact the developer.") # nolint
}

.createecdfmat <- function(scoremat, rounding, transtable, direction) {

  ecdfmat <- apply(scoremat, 2, function(x, rounding, coordvec) {
    extendedcoordvec <- rep(coordvec, round(x * rounding))
    fx <- ecdf(extendedcoordvec)(coordvec)
    return(fx)
  }, rounding, transtable$coord, simplify = TRUE)
  colnames(ecdfmat) <- gsub(direction, "", colnames(ecdfmat))
  colnames(ecdfmat) <- paste("Fx", colnames(ecdfmat), sep = "_") # nolint
  return(ecdfmat)
}

.coordandfilter <- function(str, transtable, nbrows) { # nolint

  if (isTRUE(all.equal(str, "minus"))) {
    coordvec <- seq(from = nbrows, to = 1, by = -1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable[order(coordvec), ]
    transtable <- transtable %>% dyplr::select(!dyplr::matches("plus"))
  } else {
    coordvec <- seq(from = 1, to = nbrows, by = 1)
    transtable <- cbind(transtable, coord = coordvec)
    transtable <- transtable %>% dyplr::select(!dyplr::matches("minus"))
  }
  return(transtable)
}

.computeecdf <- function(transtable, expdf, rounding, colscorevec, nbrows) { # nolint

        ## Filters the score columns according to the strand of the transcript
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
        colscorestr <- colscorevec[which(expdf$strand == str)]

        ## Create the coordinate column and select scores having the righ
        ## orientation
        direclist <- .coordandfilter(str, transtable, nbrows)
        transtable <- direclist[[1]]
        directionfill <- direclist[[2]]

        ## Building a matrix containing only the scores in the right direction
        ## for each experiment. Filling the NA values with tidyr::fill.
        scoremat <- transtable[, colscorestr]
        scoremat <- scoremat %>% tidyr::fill(contains("score"), # nolint
          .direction = directionfill)

        ## Replace the scores of transtable with the filled one
        transtable[, colscorestr] <- scoremat

        ## Retrieving the direction (fwd or rev) according to the transcript
        ## strand. It will be used to change the column names of scoremat and
        ## remove columns from transtab.
        direction <- unique(expdf[which(expdf$strand == str), "direction"])
        .checkunique(direction, "direction")
        opposedirect <- unique(expdf[which(expdf$strand != str), "direction"])
        .checkunique(opposedirect, "opposite direction") # nolint

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- .createecdfmat(scoremat, rounding, transtable, direction)

        ## Remove opposite strand from transtable and erase strand substring
        transtable <- transtable[, -grep(opposedirect, colnames(transtable))]
        colnames(transtable) <- gsub(direction, "", colnames(transtable))

        res <- cbind(transtable, ecdfmat)
        return(res)
}


#allexprsdfs=niccode_allexprsdfsvic;rounding = 10; nbcpu = nbcputrans; verbose = TRUE
genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1, # nolint
  verbose = FALSE) {

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]

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
    colnamevec <- colnames(maintable)
    colscorevec <- colnamevec[grep("_score", colnamevec)]

    ## Computing ecdf on each transcript
    if (verbose) message("\t Computing ecdf on each transcript")
    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        colscorevec, rounding, nbrows) {
        res <- .computeecdf(transtable, expdf, rounding, colscorevec, nbrows)
        return(res)
    }, expdf, colscorevec, rounding, nbrows, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(list(concatdf, nbrows))
}

!!!!!!!!!!!!!!




niccode_resecdfvic <- genesECDF(niccode_allexprsdfsvic, expdf, nbcpu = nbcputrans, verbose = TRUE)
nbwindows <- niccode_resecdfvic[[2]]
niccode_resecdfvic <- niccode_resecdfvic[[1]]

## Reading the result of ecdf that contains the column coord that is present in
## the input table of nic
viccode_resecdfvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/cugusi2023_ECDFScores_10_200.tsv"
viccode_resecdfvic <- read.table(viccode_resecdfvicpath, sep = "\t", header = TRUE)

## Adding the coordinate column
idx <- match(niccode_allexprsdfsvic[[1]]$rowid, viccode_resecdfvic$id)
if (!isTRUE(all.equal(nrow(niccode_allexprsdfsvic[[1]]), length(idx))))
    stop("Cannot retrieve coord")
niccode_allexprsdfsviccoord <- niccode_allexprsdfsvic
niccode_allexprsdfsviccoord[[1]] <- cbind(niccode_allexprsdfsviccoord[[1]], coord = viccode_resecdfvic$coord[idx])


colnames(viccode_resecdfvic) <- c("biotype.window", "chrom", "start.window",
    "end.window", "transcript", "gene", "strand.window", "window", "rowid",
    "ctrl_rep1", "ctrl_rep2", "HS_rep1", "HS_rep2", "coord", "ctrl1_score",
    "ctrl2_score", "HS1_score", "HS2_score", "Fx_ctrl1_score", "Fx_ctrl2_score",
    "Fx_HS1_score", "Fx_HS2_score")

namecolvec <- c("ctrl_rep1", "ctrl_rep2", "HS_rep1", "HS_rep2")

idxnames <- sapply(namecolvec, grep, colnames(viccode_resecdfvic))
viccode_resecdfvic <- viccode_resecdfvic[, -idxnames]
reorderednames <- c("biotype.window", "chrom", "start.window",
    "end.window", "transcript", "gene", "strand.window", "window", "rowid",
    "ctrl1_score", "ctrl2_score", "HS1_score", "HS2_score", "coord",
    "Fx_ctrl1_score", "Fx_ctrl2_score", "Fx_HS1_score", "Fx_HS2_score")
viccode_resecdfvic <- viccode_resecdfvic[, reorderednames]

if (isTRUE(all.equal(niccode_resecdfvic, viccode_resecdfvic)))
    message("table is equal after genesECDF")

