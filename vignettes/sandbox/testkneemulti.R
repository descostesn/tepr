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

dontcompvec <- showallcomp(expdf)
dontcompvec <- dontcompvec[- c(1,2,3)]

## Enter kneemulti
expthres = 0.1; nbcpu = 5; rounding = 10
dontcompare = NULL; saveobjectpath = NA; showtime = FALSE; verbose = TRUE


!!!!!!!!!!!!!! add kneemulti in vignette and showallcomp