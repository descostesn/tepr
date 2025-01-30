library("tepr")

##################
# PARAMETERS
##################

exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"

!!!!!

df <- alldf[, c()]
write.table(df, file = "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond.tsv", sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = FALSE)
!!!!!!!

##################
# MAIN
##################

## Reading and filtering on one cond several rep
expdf <- read.csv(exptabpath, header = TRUE)
expdf <- expdf[which(expdf$condition == "ctrl10"), ]

## Reading table with one cond several rep
alldf <- read.delim(finaltabpath, header = FALSE)
