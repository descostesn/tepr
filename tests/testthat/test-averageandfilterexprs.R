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
    expm <- paste0("\n No transcript was identified as expressed. You might ",
                "want to decrease the expthres parameter. Currently all genes",
                " whose expression < ", expthres, " are removed.\n")
    expect_error(averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)

    transdf$V7[which(transdf$V7 == '+')] <- "toto"
    expthres <- 0.1
    expm <- paste0("\n\t The strand name is neither \\+ or \\- in the ",
        "transcript table alldf. If you are sure to have built alldf with the ",
        "preprocessing function, contact the developer.\n")
    expect_error(averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)
})
