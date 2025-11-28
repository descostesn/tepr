## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
tmpfoldpath <- file.path(tempdir(), "tmptepr")
windsize <- 200
genomename <- "hg38"
chromtabtest <- retrievechrom(genomename, verbose = FALSE)
allchromvec <- GenomeInfoDb::seqnames(chromtabtest)
chromtabtest <- chromtabtest[allchromvec[which(allchromvec == "chr13")], ]

## Copying bedgraphs to the current directory
expdfpre <- read.csv(exptabpath)
bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
    package = "tepr"))
expdfpre$path <- bgpathvec
write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
    quote = FALSE)
exptabpath <- "exptab-preprocessing.csv"

## Necessary result to call createtablescores
allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)
blacklisthighmap(maptrackpath, blacklistpath, exptabpath, nbcputrans = 1,
    allwindowsbed, windsize, genomename = genomename, chromtab = chromtabtest,
    verbose = FALSE)

## Calling the function to test
finaltabtest <- createtablescores(tmpfold = tmpfoldpath, exptabpath,
    savefinaltable = FALSE, verbose = FALSE)


## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "finaltab.rds",
    package="tepr"))
test_that("createtablescores works properly", {
    expect_equal(finaltabtest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling createtablescores", {

    addfile <- file.path(tmpfoldpath, paste0(expdfpre[1, 1], expdfpre[1, 2],
        expdfpre[1, 3], "-chrtest.tsv"))
    write.table(data.frame(1:10), file = addfile)
    expm <- "Unequal file counts per experiment"
    expect_error(createtablescores(tmpfold = tmpfoldpath, exptabpath,
        savefinaltable = FALSE, verbose = FALSE), regexp = expm)
    file.remove(addfile)
    
    exptabtest <- rbind(expdfpre, data.frame(condition = "toto", replicate = 1,
        direction = "forward", strand = "plus", path = "toto"))
    outfile <- file.path(tempdir(), "exptabtest.csv")
    write.csv(exptabtest, file = outfile)
    expm <- "File names do not match experiment table"
    expect_error(createtablescores(tmpfold = tmpfoldpath, exptabpath = outfile,
        savefinaltable = FALSE, verbose = FALSE), regexp = expm)
})
