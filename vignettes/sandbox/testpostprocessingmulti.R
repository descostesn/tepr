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
meanctrlthres <- 0.5
meanstressthres <- 0.5
pvaltheorythres <- 0.1
aucctrlthreshigher <- -10
aucctrlthreslower <- 15
aucstressthres <- 15
attenuatedpvalksthres <- 2
outgrouppvalksthres <- 0.2
showtime <- FALSE
verbose <- TRUE


##################
# MAIN
##################

expdf <- read.csv(exptabpath, header = TRUE)
alldf <- read.delim(finaltabpath, header = FALSE)
