library("tepr")

##################
# PARAMETERS
##################

exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
# finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"
tabonecond <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond.tsv"
tabonecondonerep <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond-onerep.tsv"

##################
# MAIN
##################

## Reading and filtering on one cond several rep
expdf <- read.csv(exptabpath, header = TRUE)
expdfonecond <- expdf[which(expdf$condition == "ctrl10"), ]
expdfonerep <- expdfonecond[c(1, 2), ]

## Reading and filtering alldf
#alldf <- read.delim(finaltabpath, header = FALSE)

# df <- alldf[, c(1:9, 18:25)]
# write.table(df, file = "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond.tsv", sep = "\t",
#     quote = FALSE, row.names = FALSE, col.names = FALSE)

# dfrep <- df[, c(1:13)]
# write.table(dfrep, file = "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond-onerep.tsv", sep = "\t",
#     quote = FALSE, row.names = FALSE, col.names = FALSE)


## Reading table with one cond several rep
df <- read.delim(tabonecond, header = FALSE)

## Reading table with one cond one rep
dfrep <- read.delim(tabonecondonerep, header = FALSE)

## Parameters to enter tepr
expdf = expdfonecond; alldf = df; expthres = 0.1; nbcpu = 5; rounding = 10
dontcompare = NULL; controlcondname = "ctrl10"; stresscondname = "HS"
replaceval = NA; pval = 0.1; significant = FALSE; windsizethres = 50
countnathres = 20; meanctrlthres = 0.5; meanstressthres = 0.5
pvaltheorythres = 0.1; aucctrlthreshigher = -10; aucctrlthreslower = 15
aucstressthres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2
showtime = FALSE; verbose = TRUE
