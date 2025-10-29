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
test_that("Errors are thrown when calling tepr and teprmulti", {

    expdftest <- rbind(expdf, data.frame(condition = "test", replicate = 1,
        direction = "forward", strand = "plus",
        path = "HS_rep1_chr13.reverse.bg"))
    expm <- paste0("\n\t There are more than two conditions in your experiment",
            " table. Use teprmulti function instead.\n")
    expect_error(tepr(expdftest, transdf, expthres, verbose = FALSE),
        regexp = expm)

    expm <- paste0("\n\nThe table of values \\(alldf\\) and the table of ",
        "experiment information \\(expdf\\) do not correspond. The first four",
        " columns of expdf should be:\n\n \\-\\- condition:ctrl ctrl ctrl ctrl",
        " HS HS HS HS\n\n \\-\\- replicate: 1 1 2 2 1 1 2 2\n\n \\-\\- ",
        "direction: forward reverse forward reverse forward reverse forward ",
        "reverse\n\n \\-\\- strand: plus minus plus minus plus minus plus",
        " minus\n\n Also make sure that the bedgraph paths are correct.\n\n")
    expect_error(teprmulti(expdftest, transdf, expthres, verbose = FALSE),
        regexp = expm)
})

