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

# source("commons.R")


##################
# PARAMETERS
##################

gencodepath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
## Note: For a complete list of blacklist names see
## ah <- AnnotationHub() # nolint
## query_data <- subset(ah, preparerclass == "excluderanges") # nolint
## print(query_data) # nolint
blacklistname <- "hg38.Kundaje.GRCh38_unified_Excludable"
outputfolder <- "/g/romebioinfo/Projects/tepr/downloads"
robjoutputfold <- "/g/romebioinfo/Projects/tepr/robjsave"
## The bed file below was created and sent by Victor
maptrackpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
## Size of the window to extract values
windsize <- 200
## Table of experiments - contains the columns "name,condition,replicate,strand,path" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph.csv"
nbcpu <- 6
database_name <- "org.Hs.eg.db"


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

grepsequential <- function(valvec, gentab, invert = FALSE, verbose = FALSE) {
    invisible(sapply(valvec, function(val) {
        idx <- grep(val, gentab$V9, invert = invert)
        if (verbose)
            message(val, " - ", length(idx), " gentab - ", nrow(gentab))
        if (!isTRUE(all.equal(length(idx), 0)))
            gentab <<- gentab[idx, ]
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


bedtogr <- function(currentbed, strand = TRUE, symbol = TRUE, biotype = FALSE) {

    if (!biotype)
        grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
                ranges = IRanges::IRanges(start = currentbed[, 2],
                                      end = currentbed[, 3],
                                      names = currentbed[, 4]),
                strand = if (strand) currentbed[, 6] else "*",
                symbol = if (symbol) currentbed[, 5] else NA)
    else
        grres <- GenomicRanges::GRanges(seqnames = currentbed[, 1],
                ranges = IRanges::IRanges(start = currentbed[, 2],
                                      end = currentbed[, 3],
                                      names = currentbed[, 4]),
                strand = if (strand) currentbed[, 6] else "*",
                symbol = if (symbol) currentbed[, 5] else NA,
                biotype = currentbed[, 7])
    return(grres)
}


makewindowsbedtools <- function(expgr, binsize, biotype = FALSE) {

    ## Filtering out intervals smaller than binsize
    idxsmall <- which(GenomicRanges::width(expgr) < binsize)
    lsmall <- length(idxsmall)
    if (!isTRUE(all.equal(lsmall, 0))) {
        message("Excluding ", lsmall, "/", length(expgr), " annotations that ",
        "are too short.")
        expgr <- expgr[-idxsmall]
    }

    ## Change row names to keep the gene symbols
    if (!biotype)
        names(expgr) <- paste(names(expgr), expgr$symbol, sep = "_")
    else
        names(expgr) <- paste(names(expgr), expgr$symbol, expgr$biotype,
            sep = "_")

    ## command retrieved with HelloRanges:
    ## bedtools_makewindows("-n 200 -b stdin.bed") # nolint
    ## Note: In R, bedtools does not have the "-i srcwinnum" option
    res <- GenomicRanges::tile(expgr, n = binsize)
    res <- unlist(res, use.names = FALSE)

    ## Adding back metadata from names
    tmplist <- strsplit(names(res), "_")
    transvec <- sapply(tmplist, "[", 1)
    symbolvec <- sapply(tmplist, "[", 2)
    names(res) <- transvec
    S4Vectors::elementMetadata(res)[, "symbol"] <- symbolvec

    if (biotype) {
        biotypevec <- sapply(tmplist, "[", 3)
        S4Vectors::elementMetadata(res)[, "biotype"] <- biotypevec
    }

    ## Making names of each element of the list unique
    names(res) <- make.unique(names(res), sep = "_frame")
    return(res)
}


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


##################
# MAIN
##################

createfolder(robjoutputfold)

## Reading the information about experiments
exptab <- read.csv(exptabpath, header = TRUE)
checkexptab(exptab)

## Read gencode file
gencode <- read.delim(gencodepath, header = FALSE, skip = 5)
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (verbose) message("Combine the annotations")
protcodbed <- cbind(protcodbed, biotype = "protein-coding")
lncrnabed <- cbind(lncrnabed, biotype = "lncRNA")
allannobed <- rbind(protcodbed, lncrnabed)
allannogr <- bedtogr(allannobed, biotype = TRUE)

if (verbose) message("Make windows for all annotations")
allwindowsgr <- makewindowsbedtools(allannogr, windsize, biotype = TRUE)

if (verbose) message("Reading the black list")
blacklistgr <- createblacklist(blacklistname, outputfolder)

if (verbose) message("Reading the highly mappable ranges")
maptrack <- read.delim(maptrackpath, header = FALSE)
maptrackgr <- bedtogr(maptrack, strand = FALSE)

## Retrieving the values of the bedgraph files, removing black lists and keeping
## high mappability scores
message("Reading and filtering bedgraphs")
expnamevec <- paste0(exptab$condition, exptab$replicate, exptab$direction)
bedgraphgrlist <- retrieveandfilterfrombg(exptab, blacklistgr,
    maptrackgr, nbcpu, expnamevec)

## Retrieving values according to annotations and calculate an arithmetic
## weighted mean

# currentgr=bedgraphgrlist[[1]]
# currentstrand=exptab$strand[1]
# currentname=expnamevec[1]

## For each bedgraph
mapply(function(currentgr, currentstrand, currentname, allwindowsgr) {

    message("Overlapping ", currentname, " with annotations on strand ",
        currentstrand)
    BiocGenerics::strand(currentgr) <- currentstrand
    res <- GenomicRanges::findOverlaps(currentgr, allwindowsgr,
        ignore.strand = FALSE)

    ## Retrieving the names and frame of the mapped annotations
    idxanno <- S4Vectors::subjectHits(res)
    message("\t Retrieving transcript name and frame number")
    rownamelist <- strsplit(names(allwindowsgr)[idxanno], "_")
    transcriptvec <- sapply(rownamelist, "[", 1)
    framevec <- as.numeric(gsub("frame", "", sapply(rownamelist, "[", 2)))
    ## The frame numbering starting at 0 (empty suffix), need to transform to +1
    framevec[which(is.na(framevec))] <- 0
    framevec <- framevec + 1

    message("\t Building scoring results by transcript")
    ## Correspondance of bg score index with the transcript frame
    idxbgscorevec <- S4Vectors::queryHits(res)
    idxframedf <- data.frame(idxbgscore = idxbgscorevec, transframe = framevec,
        annoidx = idxanno)
    ## Separating the bedgraph score indexes by transcript names
    idxbgscorebytrans <- split(idxframedf, factor(transcriptvec))

    ## For each transcript, retrieve the information and the bedgraph
    ## coordinates, strand and scores, applying a weighted mean
    # tab <- idxbgscorebytrans[[1]]
    # nametrs <- names(idxbgscorebytrans)[1]
    # annogr <- allwindowsgr
    # bggr <- currentgr
    # strd <- currentstrand
    # expname <- currentname
    mapply(function(tab, nametrs, annogr, bggr, strd, expname) {

        ## Retrieving information about tables
        names(annogr) <- NULL
        annodf <- as.data.frame(annogr[tab$annoidx])
        colnames(annodf) <- paste0("trs_", colnames(annodf))
        bgdf <- as.data.frame(bggr[tab$idxbgscore])
        colnames(bgdf) <- paste0(expname, colnames(bgdf))

        ## Building the complete data.frame
        df <- do.call("cbind", list(annodf, bgdf, transcript = nametrs,
            frame = tab$transframe))

        ###########
        ## Applying a weighted mean on duplicated frames
        ###########
        dupframeidx <- which(duplicated(df$frame))
        ## For each duplicated frame
        wmeanvec <- sapply(dupframeidx, function(idxdup, df, expname) {

            ## Selecting all rows having a duplicated frame found at index idx
            allframedf <- df[which(df$frame == df$frame[idxdup]), ]
            if (isTRUE(all.equal(nrow(allframedf), 1)))
                stop("There should be more than one frame selected")

            ## Testing that the coord of the window is the same for all scores selected (this should not give an error) # nolint
            if (!isTRUE(all.equal(length(unique(allframedf[, "trs_start"])), 1)) || !isTRUE(all.equal(length(unique(allframedf[, "trs_end"])), 1))) # nolint
                stop("The size of the window is not unique for the frame rows selected, this should not happen, contact the developper.") # nolint

            ## Retrieving the coordinates and the size of the transcript
            windowstart <- allframedf[1, "trs_start"]
            windowend <- allframedf[1, "trs_end"]
            lwindow <- windowend - windowstart

            ## Retrieve the nb of overlapping nt for each score
            overntvec <- apply(allframedf, 1, function(x, expname, windowstart, windowend) {
                nt <- seq(from = x[paste0(expname, "start")], to = x[paste0(expname, "end")], by = 1)
                overnt <- length(which(nt >= windowstart & nt <= windowend))
                return(overnt)
            }, expname, windowstart, windowend)

            ## Computing weighted mean
            wmean <- weighted.mean(allframedf[, paste0(expname, "score")], overntvec)
            return(wmean)
        }, df, expname)

        ## Remove duplicated frames and replace scores by wmean

    }, idxbgscorebytrans, names(idxbgscorebytrans),
        MoreArgs = list(allwindowsgr, currentgr, currentstrand, currentname))

}, bedgraphgrlist, exptab$strand, expnamevec, MoreArgs = list(allwindowsgr),
    SIMPLIFY = FALSE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# !!!!!!!!!!!!!!!!!
# CODE TO RETRIEVE A WORKIND EXAMPLE
#
# bigtsvpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/dTAG_Cugusi_stranded_20230810.tsv" # nolint
# shdf <- read.delim(bigtsvpath, header = FALSE)
# colnames(shdf) <- c("biotype", "chr", "start", "end", "transcript", "gene", "strand", "window", "id", "name1", "score1", "name2", "score2", "name3", "score3", "name4", "score4", "name5", "score5", "name6", "score6", "name7", "score7", "name8", "score8")
# head(shdf)
# testsh <- shdf[which(shdf$gene == "ARF5"), ]
# testsh[6,]
#                biotype  chr     start       end         transcript gene strand
# 6021606 protein-coding chr7 127588491 127588507 ENST00000000233.10 ARF5      +
#         window   score1
# 6021606      6 0.440169
#
# test <- df[which(df$trs_symbol == "ARF5"), ]
# idxnodupnoz <- which(!duplicated(test$frame) & test$ctrl1fwdscore != 0)
# head(df[idxnodupnoz,])
# test <- df[which(df$frame == 6), ]
#
# !!!!!!!!!!!!!
# TEST WITH NB OF NUCLEOTIDES
# [(Nombre de NT avec coverage A)*(coverage A) + (Nombre de NT avec coverage B)*(coverage B) + … + (Nombre de NT avec coverage Z)*(coverage Z) ] / (nombre total de nucleotide pour la window) ( meme si il y a des fragments avec coverage de 0). 
#
# ## Retrieving the coordinates of the window
# windowstart <- test[1, "trs_start"]
# windowend <- test[1, "trs_end"]
# lwindow <- windowend - windowstart
#
# ## Testing that the coord of the window is the same for all scores selected (this should not give an error)
# if (!isTRUE(all.equal(length(unique(windowstart)), 1)) || !isTRUE(all.equal(length(unique(windowend)), 1)))
#     stop("The size of the window is not unique for the frame rows selected, this should not happen, contact the developper")
#
# ## Retrieve the coordinates of the scores
# nt1 <- seq(from = test[1, "ctrl1fwdstart"], to = test[1, "ctrl1fwdend"], by = 1)
# nt2 <- seq(from = test[2, "ctrl1fwdstart"], to = test[2, "ctrl1fwdend"], by = 1)
#
# ## Calculate the number of nucleotides (nt) overlapping the window
# overnt1 <- length(which(nt1 >= windowstart & nt1 <= windowend))
# overnt2 <- length(which(nt2 >= windowstart & nt2 <= windowend))
#
# ## Retrieve a percentage of overlap
# percent1 <- (100*overnt1)/lwindow
# percent2 <- (100*overnt2)/lwindow
#
# ## Retrieving scores (to be inserted in the below formula)
# score1 <- test[1, "ctrl1fwdscore"]
# score2 <- test[2, "ctrl1fwdscore"]
#
# ## Perform the weighted mean on nb of nucleotides
# ((score1*overnt1) + (score2*overnt2))/lwindow
# 0.5365867
# ## Perform the weigthed mean on nb of nucleotides with function
# weighted.mean(c(score1, score2), (overnt1, overnt2))
# 0.50305
#
# ## Perform the weighted mean on percent of nucleotides
# ((score1*percent1) + (score2*percent2))/lwindow
#  3.577244
# ## Perform the weigthed mean on nb of nucleotides with function
# weighted.mean(c(score1, score2), (percent1, percent2))
# 0.50305
#
#
# > head(testsh[,c(1:8,11)],20)
#                biotype  chr     start       end         transcript gene strand
# 6021601 protein-coding chr7 127588411 127588427 ENST00000000233.10 ARF5      +
# 6021602 protein-coding chr7 127588427 127588443 ENST00000000233.10 ARF5      +
# 6021603 protein-coding chr7 127588443 127588459 ENST00000000233.10 ARF5      +
# 6021604 protein-coding chr7 127588459 127588475 ENST00000000233.10 ARF5      +
# 6021605 protein-coding chr7 127588475 127588491 ENST00000000233.10 ARF5      +
# 6021606 protein-coding chr7 127588491 127588507 ENST00000000233.10 ARF5      +
# 6021607 protein-coding chr7 127588507 127588523 ENST00000000233.10 ARF5      +
# 6021608 protein-coding chr7 127588523 127588539 ENST00000000233.10 ARF5      +
# 6021609 protein-coding chr7 127588539 127588555 ENST00000000233.10 ARF5      +
# 6021610 protein-coding chr7 127588555 127588571 ENST00000000233.10 ARF5      +
# 6021611 protein-coding chr7 127588571 127588587 ENST00000000233.10 ARF5      +
# 6021612 protein-coding chr7 127588587 127588603 ENST00000000233.10 ARF5      +
# 6021613 protein-coding chr7 127588603 127588619 ENST00000000233.10 ARF5      +
# 6021614 protein-coding chr7 127588619 127588635 ENST00000000233.10 ARF5      +
# 6021615 protein-coding chr7 127588635 127588651 ENST00000000233.10 ARF5      +
# 6021616 protein-coding chr7 127588651 127588667 ENST00000000233.10 ARF5      +
# 6021617 protein-coding chr7 127588667 127588683 ENST00000000233.10 ARF5      +
# 6021618 protein-coding chr7 127588683 127588699 ENST00000000233.10 ARF5      +
# 6021619 protein-coding chr7 127588699 127588715 ENST00000000233.10 ARF5      +
# 6021620 protein-coding chr7 127588715 127588731 ENST00000000233.10 ARF5      +
#         window   score1
# 6021601      1 0.000000
# 6021602      2 0.000000
# 6021603      3 0.000000
# 6021604      4 0.000000
# 6021605      5 0.000000
# 6021606      6 0.440169
# 6021607      7 1.006100
# 6021608      8 1.006100
# 6021609      9 1.320506
# 6021610     10 1.949319
# 6021611     11 1.006100
# 6021612     12 1.006100
# 6021613     13 1.572031
# 6021614     14 2.012200
# 6021615     15 1.320506
# 6021616     16 1.446269
# 6021617     17 2.012200
# 6021618     18 1.446269
# 6021619     19 1.006100
# 6021620     20 1.006100
#
#    trs_seqnames trs_start   trs_end trs_width trs_strand trs_symbol
# 1          chr7 127588411 127588426        16          +       ARF5
# 2          chr7 127588427 127588442        16          +       ARF5
# 3          chr7 127588443 127588459        17          +       ARF5
# 4          chr7 127588460 127588475        16          +       ARF5
# 5          chr7 127588476 127588492        17          +       ARF5
# 6          chr7 127588493 127588508        16          +       ARF5
# 7          chr7 127588493 127588508        16          +       ARF5
# 8          chr7 127588509 127588525        17          +       ARF5
# 9          chr7 127588526 127588541        16          +       ARF5
# 10         chr7 127588542 127588558        17          +       ARF5
# 11         chr7 127588542 127588558        17          +       ARF5
# 12         chr7 127588559 127588574        16          +       ARF5
# 13         chr7 127588559 127588574        16          +       ARF5
# 14         chr7 127588575 127588590        16          +       ARF5
# 15         chr7 127588591 127588607        17          +       ARF5
# 16         chr7 127588608 127588623        16          +       ARF5
# 17         chr7 127588608 127588623        16          +       ARF5
# 18         chr7 127588624 127588640        17          +       ARF5
# 19         chr7 127588641 127588656        16          +       ARF5
# 20         chr7 127588657 127588673        17          +       ARF5
#    ctrl1fwdseqnames ctrl1fwdstart ctrl1fwdend ctrl1fwdwidth ctrl1fwdstrand
# 1              chr7     127586671   127588500          1830              +
# 2              chr7     127586671   127588500          1830              +
# 3              chr7     127586671   127588500          1830              +
# 4              chr7     127586671   127588500          1830              +
# 5              chr7     127586671   127588500          1830              +
# 6              chr7     127586671   127588500          1830              +
# 7              chr7     127588501   127588550            50              +
# 8              chr7     127588501   127588550            50              +
# 9              chr7     127588501   127588550            50              +
# 10             chr7     127588501   127588550            50              +
# 11             chr7     127588551   127588570            20              +
# 12             chr7     127588551   127588570            20              +
# 13             chr7     127588571   127588610            40              +
# 14             chr7     127588571   127588610            40              +
# 15             chr7     127588571   127588610            40              +
# 16             chr7     127588571   127588610            40              +
# 17             chr7     127588611   127588640            30              +
# 18             chr7     127588611   127588640            30              +
# 19             chr7     127588641   127588660            20              +
# 20             chr7     127588641   127588660            20              +
#    ctrl1fwdscore         transcript frame
# 1         0.0000 ENST00000000233.10     1
# 2         0.0000 ENST00000000233.10     2
# 3         0.0000 ENST00000000233.10     3
# 4         0.0000 ENST00000000233.10     4
# 5         0.0000 ENST00000000233.10     5
# 6         0.0000 ENST00000000233.10     6
# 7         1.0061 ENST00000000233.10     6
# 8         1.0061 ENST00000000233.10     7
# 9         1.0061 ENST00000000233.10     8
# 10        1.0061 ENST00000000233.10     9
# 11        2.0122 ENST00000000233.10     9
# 12        2.0122 ENST00000000233.10    10
# 13        1.0061 ENST00000000233.10    10
# 14        1.0061 ENST00000000233.10    11
# 15        1.0061 ENST00000000233.10    12
# 16        1.0061 ENST00000000233.10    13
# 17        2.0122 ENST00000000233.10    13
# 18        2.0122 ENST00000000233.10    14
# 19        1.0061 ENST00000000233.10    15
# 20        1.0061 ENST00000000233.10    16
#
# ctrl1fwdscore         transcript frame
# 1         0.0000 ENST00000000233.10     1
# 2         0.0000 ENST00000000233.10     2
# 3         0.0000 ENST00000000233.10     3
# 4         0.0000 ENST00000000233.10     4
# 5         0.0000 ENST00000000233.10     5
# 6         0.0000 ENST00000000233.10     6
# 7         1.0061 ENST00000000233.10     6
# 8         1.0061 ENST00000000233.10     7
# 9         1.0061 ENST00000000233.10     8
# 10        1.0061 ENST00000000233.10     9
# 11        2.0122 ENST00000000233.10     9
# 12        2.0122 ENST00000000233.10    10
# 13        1.0061 ENST00000000233.10    10
# 14        1.0061 ENST00000000233.10    11
# 15        1.0061 ENST00000000233.10    12
# 16        1.0061 ENST00000000233.10    13
# 17        2.0122 ENST00000000233.10    13
# 18        2.0122 ENST00000000233.10    14
# 19        1.0061 ENST00000000233.10    15
# 20        1.0061 ENST00000000233.10    16
#         window   score1
# 6021601      1 0.000000
# 6021602      2 0.000000
# 6021603      3 0.000000
# 6021604      4 0.000000
# 6021605      5 0.000000
# 6021606      6 0.440169
# 6021607      7 1.006100
# 6021608      8 1.006100
# 6021609      9 1.320506
# 6021610     10 1.949319
# 6021611     11 1.006100
# 6021612     12 1.006100
# 6021613     13 1.572031
# 6021614     14 2.012200
# 6021615     15 1.320506
# 6021616     16 1.446269
# 6021617     17 2.012200
# 6021618     18 1.446269
# 6021619     19 1.006100
# 6021620     20 1.006100
#
#
# CODE TESTED BEFORE VICTOR FORMULA (GAVE THE SAME NUMBERS)
#
# !!!!!!!!!!!!!!!
# lelement1 <- df$ctrl1fwdend[6]-df$ctrl1fwdstart[6]
#         elementcoord1 <- seq(df$ctrl1fwdstart[6], df$ctrl1fwdend[6], by = 1)
#         over1 <- which(elementcoord >= 127588491 & elementcoord <= 127588507)
#         weight1 <- (100*length(over1))/lelement1
        
# lelement2 <- df$ctrl1fwdend[7]-df$ctrl1fwdstart[7]
#         elementcoord1 <- seq(df$ctrl1fwdstart[7], df$ctrl1fwdend[7], by = 1)
#         over2 <- which(elementcoord >= 127588491 & elementcoord <= 127588507)
#         weight2 <- (100*length(over2))/lelement2
        
# !!!!!!!!!!!!!!!!!!!!!

#         lwindow <- 127588507-127588491
#         windowcoord <- seq(127588491, 127588507, by = 1)
#         idxscore1 <- which(windowcoord >= df$ctrl1fwdstart[6] & windowcoord <= df$ctrl1fwdend[6])
#         idxscore2 <- which(windowcoord >= df$ctrl1fwdstart[7] & windowcoord <= df$ctrl1fwdend[7])
#         weight1 <- (100*length(idxscore1))/lelement
#         weight2 <- (100*length(idxscore2))/lelement
#         windscore <- rep(NA, lwindow)
#         windweight <- rep(0, lwindow)
#         windscore[idxscore1] <- df$ctrl1fwdscore[6]
#         windscore[idxscore2] <- df$ctrl1fwdscore[7]
#         windweight[idxscore1] <- weight1
#         windweight[idxscore2] <- weight2

#         weighted.mean(windscore, windweight)
        # lelement1 <- df$ctrl1fwdend[6]-df$ctrl1fwdstart[6]
        # elementcoord1 <- seq(df$ctrl1fwdstart[6], df$ctrl1fwdend[6], by = 1)
        # over1 <- which(elementcoord >= df$trs_start[6] & elementcoord <= df$trs_end[6])
        # weight1 <- (100*length(over1))/lelement1
        
        # lelement2 <- df$ctrl1fwdend[7]-df$ctrl1fwdstart[7]
        # elementcoord1 <- seq(df$ctrl1fwdstart[7], df$ctrl1fwdend[7], by = 1)
        # over2 <- which(elementcoord >= df$trs_start[6] & elementcoord <= df$trs_end[6])
        # weight2 <- (100*length(over2))/lelement2
        

        # lwindow <- df$trs_end[6]-df$trs_start[6]
        # windowcoord <- seq(df$trs_start[6], df$trs_end[6], by = 1)
        # idxscore1 <- which(windowcoord >= df$ctrl1fwdstart[6] & windowcoord <= df$ctrl1fwdend[6])
        # idxscore2 <- which(windowcoord >= df$ctrl1fwdstart[7] & windowcoord <= df$ctrl1fwdend[7])
        # windscore <- rep(NA, lwindow)
        # windweight <- rep(0, lwindow)
        # windscore[idxscore1] <- df$ctrl1fwdscore[6]
        # windscore[idxscore2] <- df$ctrl1fwdscore[7]
        # windweight[idxscore1] <- weight1
        # windweight[idxscore2] <- weight2

        # weighted.mean(windscore, windweight)
        # !!!!!!!!!!!!!!!!!!!



## Saving objects to check conformity with bash results
saveRDS(protcodbed, file = file.path(robjoutputfold, "protcodbed.rds"))
saveRDS(protcodgr, file = file.path(robjoutputfold, "protcodgr.rds"))
saveRDS(lncrnabed, file = file.path(robjoutputfold, "lncrnabed.rds"))
saveRDS(lncrnagr, file = file.path(robjoutputfold, "lncrnagr.rds"))
# protcodbed <- readRDS(file.path(robjoutputfold, "protcodbed.rds"))
# protcodgr <- readRDS(file.path(robjoutputfold, "protcodgr.rds"))
# lncrnabed <- readRDS(file.path(robjoutputfold, "lncrnabed.rds"))
# lncrnagr <- readRDS(file.path(robjoutputfold, "lncrnagr.rds"))


## Exclude blacklist
blacklistgr <- createblacklist(blacklistname, outputfolder)
protcodnoblackgr <- excludeorkeepgrlist(protcodgr, blacklistgr)
lncrnanoblackgr <- excludeorkeepgrlist(lncrnagr, blacklistgr)

## Check excluded intervals using blacklist
checkremoval(protcodgr, protcodnoblackgr, "proteincoding", "blacklist",
    blacklistgr, removeopt = TRUE)
checkremoval(lncrnagr, lncrnanoblackgr, "lncrna", "blacklist",
    blacklistgr, removeopt = TRUE)

## Exclude low mappability
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES - the mappability track
## used has only 1 as mapping scores. See parameters.
maptrack <- read.delim(maptrackpath, header = FALSE)
maptrackgr <- bedtogr(maptrack)


!!!!!!!!
keepgrlist <- function(datagr, maptrackgr) {
    res <- GenomicRanges::findOverlaps(datagr, maptrackgr)
    idxtokeep <- unique(S4Vectors::queryHits(res))
    datagr <- datagr[idxtokeep, ]
    return(datagr)
}
protcodnoblackmaponlygr <- keepgrlist(protcodnoblackgr, maptrackgr)
lncrnanoblackmaponlygr <- keepgrlist(lncrnanoblackgr, maptrackgr)

!!!!!!!!!!!

protcodnoblacknomapgr <- excludeorkeepgrlist(protcodnoblackgr, maptrackgr,
    removefrom = FALSE)
lncrnanoblacknomapgr <- excludeorkeepgrlist(lncrnanoblackgr, maptrackgr,
    removefrom = FALSE)

## Check excluded intervals because of low mappability for protein coding
checkremoval(protcodnoblackgr, protcodnoblacknomapgr, "proteincoding",
    "maptrack", maptrackgr, removeopt = FALSE)
checkremoval(lncrnanoblackgr, lncrnanoblacknomapgr, "lncrna", "maptrack",
    maptrackgr, removeopt = FALSE)

## Make windows of windsize for each annotation
## WARNING: CANNOT FIND EXACTLY THE SAME NUMBER OF LINES

lncrnawindows <- makewindowsbedtools(lncrnanoblacknomapgr, windsize)

## Retrieving values from bigwig files
protcoddf <- buildscoreforintervals(protcodwindows, exptab, "protein_coding",
    nbcpu, database_name)
lncrnadf <- buildscoreforintervals(lncrnawindows, exptab, "lncrna", nbcpu,
    database_name)
alldf <- rbind(protcoddf, lncrnadf)
saveRDS(alldf, file = file.path(robjoutputfold, "alldffrompreprocessing.rds"))
#alldf <- readRDS("/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds") # nolint
