library("tepr")

##################
# PARAMETERS
##################

exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"
expthres <- 0.1
nbcpu <- 5
rounding <- 10
dontcompare <- NULL
controlcondname <- "ctrl"
stresscondname <- "HS"
replaceval <- NA
pval <- 0.1
significant <- FALSE
windsizethres <- 50
countnathres <- 20
meancond1thres <- 0.5
meancond2thres <- 0.5
pvaltheorythres <- 0.1
auccond1threshigher <- -10
auccond1threslower <- 15
auccond2thres <- 15
attenuatedpvalksthres <- 2
outgrouppvalksthres <- 0.2
showtime <- FALSE
verbose <- TRUE


##################
# MAIN
##################

expdf <- read.csv(exptabpath, header = TRUE)
alldf <- read.delim(finaltabpath, header = FALSE)
