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
ecdf <- genesECDF(avfilt, verbose = FALSE)
resecdf <- ecdf[[1]]
nbwindows <- ecdf[[2]]

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "meandiff.rds",
    package="tepr"))
meandifftest <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
test_that("meandifference works properly", {
    expect_equal(meandifftest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling meandifference", {

    expdftest <- expdf
    expdftest$condition[which(expdf$condition == "ctrl")] <- "toto"
    expm <- "Condition not found"
    expect_error(meandifference(resecdf, expdftest, nbwindows, verbose = FALSE),
        regexp = expm)
    
    resecdftest <- resecdf
    idxfx <- grep("Fx", colnames(resecdftest))
    colvecfx <- colnames(resecdftest)[idxfx]
    newcolvec <- gsub("Fx", "fx", colvecfx)
    colnames(resecdftest)[idxfx] <- newcolvec
    expm <- "Missing columns in meandifference"
    expect_error(meandifference(resecdftest, expdf, nbwindows, verbose = FALSE),
        regexp = expm)
})
