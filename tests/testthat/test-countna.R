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
             expect_identical(countnatest, expectedobj)
         })

## ----- Checking errors ----- ##
test_that("Errors are thrown when calling countna", {})

allexprsdfs, expdf, nbcpu = 1, showtime = FALSE,
  verbose = TRUE