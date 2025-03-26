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
    expect_equal(allwindowsbed, expectedobj)
})
