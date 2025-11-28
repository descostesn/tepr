## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilt <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
ecdf <- genesECDF(avfilt, verbose = FALSE)
resecdf <- ecdf[[1]]
nbwindows <- ecdf[[2]]
resmeandiff <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "attenuation.rds",
    package="tepr"))
resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
        resmeandiff, verbose = FALSE)
test_that("attenuation works properly", {
    expect_equal(resatt, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling attenuation", {

    kneetest <- resknee
    kneetest[, "knee_AUC_ctrl"] <- 0
    expm <- "Internal error retrieving idxup"
    expect_error(attenuation(resauc, kneetest, rescountna, bytranslistmean,
        expdf, resmeandiff, verbose = FALSE), regexp = expm)
    
    kneetest <- resknee
    kneetest[, "knee_AUC_ctrl"] <- 400
    expm <- "Internal error retrieving idxdown"
    expect_error(attenuation(resauc, kneetest, rescountna, bytranslistmean,
        expdf, resmeandiff, verbose = FALSE), regexp = expm)
})
