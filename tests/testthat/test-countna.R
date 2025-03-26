## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating averageandfilterexprs to call countNA
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilttest <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)

## ---- Comparing to expected object ---- ##
expectedobj <- readRDS(system.file("extdata", "countna.rds",
    package="tepr"))
countnatest <- countna(avfilttest, expdf, nbcpu = 1, verbose = FALSE)
test_that("countna works properly", {
    expect_equal(countnatest, expectedobj)
})

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling countna", {

    avfilt <- avfilttest
    avfilt[[1]]$strand[which(avfilt[[1]]$strand == '+')] <- "toto"
    expm <- paste0("\n\t In .computeecdf or countna, strand is neither plus or",
            " minus in the table returned by the function ",
            "averageandfilterexprs. This should not happen. Contact the ",
            "developer.\n")
    expect_error(countna(avfilt, expdf, verbose = FALSE), regexp = expm)

    avfilt <- avfilttest
    avfilt[[1]][which(avfilt[[1]]$gene == "AP5S1")[c(1,2)],
        "ctrl_rep1.plus_score"] <- NA
    avfilt[[1]][which(avfilt[[1]]$gene == "AP5S1")[c(1,2)],
        "ctrl_rep2.plus_score"] <- NA
    expm <- paste0("\n\t Number of NA is different between conditions for ",
              "AP5S1: 10 - 8",
              ". This should not happen. Contact the developer.\n")
    expect_error(countna(avfilt, expdf, verbose = FALSE), regexp = expm)
})
