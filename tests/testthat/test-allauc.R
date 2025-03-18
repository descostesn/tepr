## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilt <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
countna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
resecdf <- ecdf[[1]]
nbwindows <- ecdf[[2]]
meandiff <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
bytranslistmean <- split(meandiff, factor(meandiff$transcript))

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "allauc.rds",
    package="tepr"))
allauctest <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
test_that("allauc works properly", {
             expect_identical(allauctest, expectedobj)
         })
