## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package = "tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
  package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"

## Copying bedgraphs to the current directory
expdfpre <- read.csv(exptabpath)
bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
    package = "tepr"))
expdfpre$path <- bgpathvec
write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
    quote = FALSE)
exptabpath <- "exptab-preprocessing.csv"

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "finaltab.rds",
    package="tepr"))
finaltabtest <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = genomename, verbose = FALSE)
test_that("preprocessing works properly", {
    expect_equal(finaltabtest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling preprocessing", {

    expm <- "\n\t Either the genome name or chromtab should be provided.\n"
    expect_error(preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = NA), regexp = expm)

    rdsfile <- file.path(tempdir(), "finaltable.rds")
    saveRDS(1, file = rdsfile)
    expm <- paste0("\n\t The final table already exists, set reload = FALSE to",
            " create it again.\n")
    expect_error(preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = "hg38", reload = TRUE), regexp = expm)
    file.remove(rdsfile)

    expm <- paste0("\n Chromtab should be a Seqinfo object. Use ",
        "rtracklayer::SeqinfoForUCSCGenome\\(genomename\\).\n")
    expect_error(preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = "hg38", chromtab = 2), regexp = expm)
    
    chromtabtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
    expm <- paste0("\n Non-canonical chromosomes found in chromtab. If you",
                    " are sure you want to proceed set forcechrom = TRUE.\n\n")
    expect_error(suppressWarnings(preprocessing(exptabpath, gencodepath,
        windsize, maptrackpath, blacklistpath, genomename = "hg38",
        chromtab = chromtabtest)), regexp = expm)    
})
