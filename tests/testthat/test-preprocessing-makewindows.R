## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
windsize <- 200

## Necessary result to call makewindows
allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)

## Calling the function to test
allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "allwindowsbed.rds",
    package="tepr"))
test_that("makewindows works properly", {

    ## Darwin indicates a macos. Because of floating system, object is
    ## equal but not identical on this OS.
    if (!isTRUE(all.equal(Sys.info()[['sysname']], "Darwin"))) {
        expect_identical(allwindowsbed, expectedobj)
    } else {
        expect_equal(allwindowsbed, expectedobj)
    }
})
