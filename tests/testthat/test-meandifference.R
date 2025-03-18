## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
countnatest <- countna(avfilttest, expdf, nbcpu = 1, verbose = FALSE)
ecdftest <- genesECDF(avfilttest, expdf, verbose = FALSE)
resecdf <- ecdftest[[1]]
nbwindows <- ecdftest[[2]]

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "meandiff.rds",
    package="tepr"))
meandifftest <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
test_that("meandifference works properly", {
             expect_identical(meandifftest, expectedobj)
         })
