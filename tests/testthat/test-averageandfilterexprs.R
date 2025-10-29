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

    expthres <- 0.1
    expdfshift <- expdf[c(6, 5, 3, 8, 4, 1, 7, 2),]
    expm <- paste0("\n\nThe table of values \\(alldf\\) and the table of ",
        "experiment information \\(expdf\\) do not correspond. The first four",
        " columns of expdf should be:\n\n \\-\\- condition:ctrl ctrl ctrl ctrl",
        " HS HS HS HS\n\n \\-\\- replicate: 1 1 2 2 1 1 2 2\n\n \\-\\- ",
        "direction: forward reverse forward reverse forward reverse forward ",
        "reverse\n\n \\-\\- strand: plus minus plus minus plus minus plus",
        " minus\n\n Also make sure that the bedgraph paths are correct.\n\n")
    expect_error(averageandfilterexprs(expdfshift, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)

    transdf$V7[which(transdf$V7 == '+')] <- "toto"
    expm <- paste0("\n\t The strand name is neither \\+ or \\- in the ",
        "transcript table alldf. If you are sure to have built alldf with the ",
        "preprocessing function, contact the developer.\n")
    expect_error(averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE), regexp = expm)
})
