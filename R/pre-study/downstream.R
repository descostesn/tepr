####################################
# This script goes through documentation/explore.R and homogenizes it with
# preprocessing.R
#
# Descostes - R-4.4.1 - July 2024
####################################



##################
# PARAMETERS
##################

alldfpath <- "/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv"
expthres <- 0.1
nbcpu <- 6


##################
# MAIN
##################

## Reading alldf and info tab
alldf <- readRDS(alldfpath)
exptab <- read.csv(exptabpath, header = TRUE)

## Filtering out non expressed transcripts:
## 1) for each column, calculate the average expression per transcript (over each frame) # nolint
## 2) For each column, remove a line if it contains value < expthres
idxscores <- sapply(exptab$name, grep, colnames(alldf))
if (isTRUE(all.equal(length(idxscores), 0)))
    stop("The scores were not retrieved in the data.frame") # nolint

# idxremove <- unique(unlist(apply(alldf[, idxscores], 2,
#     function(currentcol, thres) { return(which(currentcol < thres)) },
#     expthres, simplify = FALSE)))
