########################
## This script aims at reproducing the figure 2d of the article using the
## object completedf generated with the script downstream.R
##
## Descostes - R-4.4.1 - sept 2024
########################


##################
# PARAMETERS
##################


completedfpath <- "/g/romebioinfo/tmp/downstream/completedf.rds"
dfmeandiffpath <- "/g/romebioinfo/tmp/downstream/dfmeandiff.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
colvec <- c("#90AFBB", "#FF9A04", "#10AFBB", "#FC4E07")
genename <- "EGFR"


##################
# MAIN
##################

## Reading inputs
completedf <- readRDS(completedfpath)
dfmeandiff <- readRDS(dfmeandiffpath)
expdf <- read.csv(exptabpath, header = TRUE)

plotecdf <- function(dfmeandiff, completedf, genename, verbose = TRUE) {

    ## Retrieving rows concerning the gene of interest
    if (verbose) message("\t Retrieving rows concerning the gene of interest")
    df <- dfmeandiff[which(dfmeandiff$gene == genename), ]

    ## Retrieving gene info
    geneinfo <- completedf[which(completedf$gene == genename), ]

    ## Computing the window size factor
    df100 <- df[which(df$coord == 100), ]
    windsizefact <- (df100$end - df100$start) / 1000
}
