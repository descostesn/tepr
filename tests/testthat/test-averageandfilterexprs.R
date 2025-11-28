## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Read input tables
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "averageandfilterexprs.rds",
    package="tepr"))
avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
test_that("preprocessing works properly", {
    expect_equal(avfilttest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling averageandfilterexprs", {

    expthres <- 10000
    expm <- "No expressed transcripts found"
    expect_error(averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)

    expthres <- 0.1
    expdfshift <- expdf[c(6, 5, 3, 8, 4, 1, 7, 2),]
    expm <- "Experiment table mismatch"
    expect_error(averageandfilterexprs(expdfshift, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)

    transdf$V7[which(transdf$V7 == '+')] <- "toto"
    expm <- "Invalid strand value"
    expect_error(averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)
})
