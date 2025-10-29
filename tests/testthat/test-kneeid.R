## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilt <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
resecdf <- ecdf[[1]]
nbwindows <- ecdf[[2]]
meandiff <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
bytranslistmean <- split(meandiff, factor(meandiff$transcript))

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "kneeid.rds",
    package="tepr"))
kneeidtest <- kneeid(bytranslistmean, expdf, verbose = FALSE)
test_that("kneeid works properly", {
    expect_equal(kneeidtest, expectedobj)
})
