## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
expdf <- read.csv(exppath)

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling checkexptab", {

    exptest <- expdf
    colnames(exptest) <- c("cond", "replicate", "direction", "strand", "path")
    expm <- paste0("The experiment table should have the columns: 'condition',",
        " 'replicate', 'direction', 'strand', 'path'.")
    expect_error(checkexptab(exptest), regexp = expm)
})

