## Parameters
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.0000000000000000000000001

## Read input tables
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)

reslist <- averageandfilterexprs(expdf, transdf, expthres, showtime = FALSE, # nolint
    verbose = TRUE)
