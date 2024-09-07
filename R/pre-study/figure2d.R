########################
## This script aims at reproducing the figure 2d of the article using the
## object completedf generated with the script downstream.R
##
## Descostes - R-4.4.1 - sept 2024
########################

library("ggplot2")
library("tidyr")
library("tidyselect")


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

!!!!!!!!!!!!

.subtext <- function(condvec, geneinfo, digits) {

    subtext <- sapply(condvec, function(cond, geneinfo, digits) {

            ksname <- paste0("adjFDR_pvalaucks_", cond)
            ksval <- round(geneinfo[, ksname], digits)
            aucval <- round(geneinfo[, paste0("auc_", cond)], digits)

            if (ksval < 0.01) {
                kneeaucval <- geneinfo[, paste0("knee_AUC_", cond)]
                kneeval <- round(kneeaucval * geneinfo$windsize / 1000, digits)
            } else {
                kneeval <- "NA"
            }
            return(paste0(cond, ": AUC = ", aucval, ", KS = ", ksval,
                ", Knee (kb) = ", kneeval))
        }, geneinfo, digits)

    return(subtext)
}

.valcolbuild <- function(condvec, repvec) {
    return(as.vector(sapply(condvec,
        function(cond) { sapply(repvec, function(rep) {
            paste0(cond, rep, "score") })})))
}

plotecdf <- function(dfmeandiff, completedf, genename, digits = 2,
    verbose = TRUE) {

    ## Retrieving rows concerning the gene of interest
    if (verbose) message("\t Retrieving rows concerning the gene of interest")
    df <- dfmeandiff[which(dfmeandiff$gene == genename), ]
    geneinfo <- completedf[which(completedf$gene == genename), ]

    if (verbose) message("\t Gathering statistics about each condition")
    ## Computing the window size factor
    df100 <- df[which(df$coord == 100), ]
    windsizefact <- (df100$end - df100$start) / 1000
    ## Retrieving auc, ks, and knee
    condvec <- unique(expdf$condition)
    subtext <- .subtext(condvec, geneinfo, digits)

    ## Building data.frame for plotting with fx and value
    repvec <- unique(expdf$replicate)
    colnamedfvec <- colnames(df)
    fxcolvec <- colnamedfvec[grep("^Fx_", colnamedfvec)]
    valcolvec <- .valcolbuild(condvec, repvec)
    ## Apply pivot
    dflongfx <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(fxcolvec), names_to = "conditions",
        values_to = "Fx") %>%
        dplyr::mutate(conditions = gsub("Fx_|_score", "", conditions)) # nolint
    dflongval <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(valcolvec), names_to = "conditions",
        values_to = "value") %>%
        dplyr::mutate(conditions = gsub("value_|_score", "", conditions)) # nolint

}
