## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")

## Call the function to test
allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "allannobed.rds",
    package="tepr"))
test_that("retrieveanno works properly", {
    expect_equal(allannobed, expectedobj)
})
