## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "tepr.rds",
    package="tepr"))
restepr <- tepr(expdf, transdf, expthres, verbose = FALSE)
test_that("tepr works properly", {
    expect_equal(restepr, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling tepr", {

    expdftest <- rbind(expdf, data.frame(condition = "test", replicate = 1,
        direction = "forward", strand = "plus",
        path = "HS_rep1_chr13.reverse.bg"))
    expm <- paste0("\n\t There are more than two conditions in your experiment",
            " table. Use teprmulti function instead.\n")
    expect_error(tepr(expdftest, transdf, expthres, verbose = FALSE),
        regexp = expm)
})
