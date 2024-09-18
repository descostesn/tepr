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
rm(allbgnic)
gc()

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
    nbcpubg, allwindowsbed, expnamevec, windsize, verbose = TRUE) {

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
    #  currentpath=exptab$path[1];currentname=expnamevec[1];
    #     currentstrand=exptab$strand[1]
    bedgraphlist <- parallel::mcmapply(function(currentpath, currentname,
        currentstrand, allwindtib, blacklisttib, maptracktib, nbcpuchrom,
        windsize, verbose) {

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
        suppressWarnings(resanno <- valr::bed_intersect(valtib, allwindstrand,
            suffix = c("", ".anno")))

        ## Removing black list
        if (verbose) message("\t Keeping scores not on black list")
        resblack <- valr::bed_intersect(resanno, blacklisttib, invert = TRUE)

        ## Processing by chromosomes because of size limits, the mappability
        ## track has too many rows
        if (verbose) message("\t Keeping scores on high mappability track")
        chromvec <- as.data.frame(unique(maptracktib["chrom"]))[, 1]
        resmaplist <- lapply(chromvec, function(currentchrom, allwindtib) {

            if (verbose) message("\t\t over ", currentchrom)
            ## Keeping scores on high mappability track
            resmap <-  valr::bed_intersect(resblack,
                maptracktib  %>% dplyr::filter(chrom == currentchrom), # nolint
                suffix = c("", "maphigh"))

            ## Removing mapping columns and duplicates
            message("\t\t\t Removing mapping columns and duplicates")
            resmap <- resmap[,-grep("maphigh|.overlap|.source",
                colnames(resmap))]
            resmap <- resmap %>% dplyr::distinct(chrom, start, end,
                start.anno, end.anno, .keep_all = TRUE)
            invisible(gc())

            ## Processing data per transcript
            message("\t\t\t Building scoring results by transcript")
            bgscorebytrans <- split(resmap, factor(resmap$transcript.anno))

            #currenttrans=bgscorebytrans[[1]]
            lapply(bgscorebytrans, function(currenttrans, windsize, resmap) {
                ## Setting missing frames to NA
               idx <- match(seq_len(windsize), unique(currenttrans$window.anno))
               idxna <- which(is.na(idx))
               if (!isTRUE(all.equal(length(idxna), 0))) {
                    message("\t\t\t Integrating missing scores")
                    missinglist <- lapply(idxna, function(currentna, resmap,
                        currenttrans) {
                            misschrom <- currenttrans$chrom[1]
                            misstrans <- currenttrans$transcript.anno[1]
                            missgene <- currenttrans$gene.anno[1]

                            idxmiss <- which(resmap$chrom == misschrom &
                                resmap$transcript.anno == misstrans &
                                resmap$gene.anno == missgene)
                            
                            !!!!!!!!! RETRIEVE CURRENTNA - 1 EXCEPT IF CURRENTNA IS 1
                            !!!!!!!!!!! THEN MODIFY CURRENTTRANS WITH <<-
                            !!!!!!!!!! CHANGE LOOP TO INVISIBLE
                            if (!isTRUE(all.equal(length(idxmiss), 1)))
                                stop("idxmiss should be unique in retrieveandfilterfrombg. Contact the developper.")
                            
                            resrow <- resmap[idxmiss, ]
                            resrow["score"] <- NA

                              chrom     start       end width strand score   biotype.anno start.anno
1  chr7 127586671 127588500  1830      *     0 protein-coding  127588411
2  chr7 127586671 127588500  1830      *     0 protein-coding  127588427
   end.anno    transcript.anno gene.anno strand.anno window.anno coord.anno
1 127588427 ENST00000000233.10      ARF5           +           1          1
2 127588443 ENST00000000233.10      ARF5           +           2          2

                            chrom  start    end width strand  score biotype.anno start.anno end.anno
1  chr7 149501 149640   140      * 0.0000       lncRNA     149597   149626
2  chr7 149501 149640   140      * 0.0000       lncRNA     149626   149655
3  chr7 149641 149670    30      * 1.0061       lncRNA     149626   149655
4  chr7 149641 149670    30      * 1.0061       lncRNA     149655   149684
5  chr7 149671 149730    60      * 2.0122       lncRNA     149655   149684
6  chr7 149671 149730    60      * 2.0122       lncRNA     149684   149713
    transcript.anno gene.anno strand.anno window.anno coord.anno
1 ENST00000484550.1 LINC03014           +           1          1
2 ENST00000484550.1 LINC03014           +           2          2
3 ENST00000484550.1 LINC03014           +           2          2
4 ENST00000484550.1 LINC03014           +           3          3
5 ENST00000484550.1 LINC03014           +           3          3
6 ENST00000484550.1 LINC03014           +           4          4
               }
                
            }, windsize, resmap)

            !!
            
            return(resmap)
            }, allwindtib)

        return(resmap)

    }, exptab$path, expnamevec, exptab$strand, MoreArgs = list(allwindtib,
        blacklisttib, maptracktib, windsize, verbose), SIMPLIFY = FALSE,
        mc.cores = nbcpubg)

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
