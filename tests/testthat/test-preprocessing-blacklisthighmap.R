## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"
chromtabtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
allchromvec <- GenomeInfoDb::seqnames(chromtabtest)


## Copying bedgraphs to the current directory
expdfpre <- read.csv(exptabpath)
bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
    package = "tepr"))
expdfpre$path <- bgpathvec
write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
    quote = FALSE)
exptabpath <- "exptab-preprocessing.csv"

## Necessary result to call blacklisthighmap
allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling blacklisthighmap", {

    expm <- "\n\t Either the genome name or chromtab should be provided.\n"
    expect_error(blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
    nbcputrans = 1, allwindowsbed, windsize), regexp = expm)

    expm <- paste0("\n Chromtab should be a Seqinfo object. Use ",
        "rtracklayer::SeqinfoForUCSCGenome\\(genomename\\).\n")
    expect_error(blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
    nbcputrans = 1, allwindowsbed, windsize, chromtab = 2), regexp = expm)

    expm <- paste0("\n Non-canonical chromosomes found in chromtab. If you",
                    " are sure you want to proceed set forcechrom = TRUE.\n\n")
    expect_error(blacklisthighmap(maptrackpath, blacklistpath,
        exptabpath, nbcputrans = 1, allwindowsbed, windsize,
        chromtab = chromtabtest, verbose = FALSE), regexp = expm)
    
    expm <- paste0("\n\t The genome toto was not found with the function",
        " rtracklayer::SeqinfoForUCSCGenome. Check the spelling or verify",
        " if the genome is available on UCSC. The connection to UCSC can ",
        "also have some hickup. You can callagain the function using the ",
        "chromtab parameter: chromtab <- rtracklayer::SeqinfoForUCSCGenome\\(",
        "genomename\\).\n")
    expect_error(blacklisthighmap(maptrackpath, blacklistpath,
        exptabpath, nbcputrans = 1, allwindowsbed, windsize,
        genomename = "toto", verbose = FALSE), regexp = expm)
})
