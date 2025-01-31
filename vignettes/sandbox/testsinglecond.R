library("tepr")

##################
# PARAMETERS
##################

## DRB data
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"
tabonecond <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond.tsv"
tabonecondonerep <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq-onecond-onerep.tsv"

## Cugusi data
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph.csv" # nolint
dTAG_Cugusi_stranded_20230810.tsv
finaltabpath <- "/g/romebioinfo/tmp/preprocessing/objects-tsv-15cpus/cugusi.tsv"
tabonecond <- "/g/romebioinfo/tmp/preprocessing/objects-tsv-15cpus/cugusi-onecond.tsv"
tabonecondonerep <- "/g/romebioinfo/tmp/preprocessing/objects-tsv-15cpus/cugusi-onecond-onerep.tsv"



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
# expdf = expdfonerep; alldf = dfrep; expthres = 0.1; nbcpu = 5; rounding = 10
# dontcompare = NULL; controlcondname = "ctrl10"; stresscondname = "HS"
# replaceval = NA; pval = 0.1; significant = FALSE; windsizethres = 50
# countnathres = 20; meanctrlthres = 0.5; meanstressthres = 0.5
# pvaltheorythres = 0.1; aucctrlthreshigher = -10; aucctrlthreslower = 15
# aucstressthres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2
# showtime = FALSE; verbose = TRUE

# load tepr on one cond
rescond <- tepr(expdf = expdfonecond, alldf = df, expthres = 0.1, nbcpu = 5,
    controlcondname = "ctrl10", showtime = TRUE, verbose = TRUE)

# load tepr on one cond one rep
resrep <- tepr(expdf = expdfonerep, alldf = dfrep, expthres = 0.1, nbcpu = 5,
    controlcondname = "ctrl10", showtime = TRUE, verbose = TRUE)

## Testing plotecdf on one cond
dfmeandiff = rescond[[1]]; unigroupdf = rescond[[2]]; expdf = expdfonecond
genename = "COQ9"; colvec = c("#90AFBB", "#10AFBB")
outfold = "."; digits = 2; middlewind = 100; pval = 0.5; plot = TRUE
formatname = "pdf"; verbose = TRUE
