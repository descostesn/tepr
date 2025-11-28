## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
expdf <- read.csv(exppath)

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling checkexptab", {

    exptest <- expdf
    colnames(exptest) <- c("cond", "replicate", "direction", "strand", "path")
    expm <- "Missing columns in experiment table"
    expect_error(checkexptab(exptest), regexp = expm)

    exptest <- expdf
    exptest[1, 3] <- "toto"
    expm <- "Invalid 'direction' values"
    expect_error(checkexptab(exptest), regexp = expm)

    exptest <- expdf
    exptest[1, 4] <- "toto"
    expm <- "Invalid 'strand' values"
    expect_error(checkexptab(exptest), regexp = expm)

    exptest <- expdf
    exptest[1, 1] <- "toto_wt"
    expm <- "Invalid condition names"
    expect_error(checkexptab(exptest), regexp = expm)
})

