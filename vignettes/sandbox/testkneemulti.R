library("tepr")

##################
# PARAMETERS
##################

## DRB data
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"


##################
# MAIN
##################


########################
## DRB
########################

## Reading and filtering alldf
alldf <- read.delim(finaltabpath, header = FALSE) # nolint
expdf <- read.csv(exptabpath, header = TRUE)

## Enter teprmulti

dontcompare <- c("ctrl0_vs_ctrl10","ctrl0_vs_ctrl20","ctrl0_vs_ctrl30","ctrl0_vs_ctrl40")



!!!!!!!!!! code utils show all comparisons
!!!!!!!!!!!!!! add kneemulti in vignette and showallcomp