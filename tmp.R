library("ggplot2")
library("tidyr")
library("tidyselect")

## /g/romebioinfo/tmp/comparewithscratch-plotting


##################
# PARAMETERS
##################

unigroupdfpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_unigroupdf.rds" # nolint
dfmeandiffpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_dfmeandiffvic.rds" # nolint
expdfpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
colvec <- c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
outfold <- "/g/romebioinfo/tmp/comparewithscratch-plotting"


##################
#FUNCTIONS
##################











##################
# MAIN
##################

## Reading objects and files
dfmeandiff <- readRDS(dfmeandiffpath)
unigroupdf <- readRDS(unigroupdfpath)
expdf <- read.csv(expdfpath, header = TRUE)


####
#### plotecdf
####

plotecdf(dfmeandiff, unigroupdf, expdf, "EGFR", colvec, outfold, plot = TRUE)
plotecdf(dfmeandiff, unigroupdf, expdf, "MARCHF6", colvec, outfold, plot = TRUE)


####
#### plotauc
####

genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")
plotauc(unigroupdf, genevec, plot = TRUE)
plotauc(unigroupdf, legendpos = "none", subtitle = "Genes selected for Unibind",
    maintitle = "AUC Control vs HS", plot = TRUE, plottype = "groups")


####
#### metagene
####

plotmetagenes(unigroupdf, "attenuation", plot = TRUE)
plotmetagenes(unigroupdf, "outgroup", plot = TRUE)
plotmetagenes(unigroupdf, "universe", plot = TRUE)
plotmetagenes(unigroupdf, "all", plot = TRUE)


####
#### histogram
####

plothistoknee(unigroupdf, plottype = "percent", plot = TRUE)
plothistoknee(unigroupdf, plottype = "kb", plot = TRUE)
