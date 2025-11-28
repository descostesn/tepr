## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"
chromtabtest <- retrievechrom(genomename, verbose = FALSE, filterchrom = FALSE)
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

    expm <- "Missing genome information"
    expect_error(blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
    nbcputrans = 1, allwindowsbed, windsize), regexp = expm)

    expm <- "Invalid chromtab type"
    expect_error(blacklisthighmap(maptrackpath, blacklistpath, exptabpath,
    nbcputrans = 1, allwindowsbed, windsize, chromtab = 2), regexp = expm)

    expm <- "Non-canonical chromosomes in chromtab"
    expect_error(blacklisthighmap(maptrackpath, blacklistpath,
        exptabpath, nbcputrans = 1, allwindowsbed, windsize,
        chromtab = chromtabtest, verbose = FALSE), regexp = expm)
    
    expm <- "Genome not found"
    expect_error(blacklisthighmap(maptrackpath, blacklistpath,
        exptabpath, nbcputrans = 1, allwindowsbed, windsize,
        genomename = "toto", verbose = FALSE), regexp = expm)
})
