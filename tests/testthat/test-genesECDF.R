## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating averageandfilterexprs and countNA to call genesECDF
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
countnatest <- countna(avfilttest, expdf, nbcpu = 1, verbose = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "genesecdf.rds",
    package="tepr"))
ecdftest <- genesECDF(avfilttest, verbose = FALSE)
test_that("genesECDF works properly", {
    expect_equal(ecdftest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling genesECDF", {

    avfilt <- avfilttest
    avfilt[[1]]$strand[which(avfilt[[1]]$strand == '+')] <- "toto"
    expm <- "Invalid strand value"
    expect_error(genesECDF(avfilt, verbose = FALSE), regexp = expm)
})
