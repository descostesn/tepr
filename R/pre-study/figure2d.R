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
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
colvec <- c("#90AFBB", "#FF9A04", "#10AFBB", "#FC4E07")
genename <- "EGFR"


##################
# MAIN
##################

## Reading inputs
completedf <- readRDS(completedfpath)
expdf <- read.csv(exptabpath, header = TRUE)

