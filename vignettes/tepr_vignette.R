## ----quickstart, results = 'hide', message = FALSE, warning = FALSE-----------

library(tepr)

#########
# Pre-processing - takes ~ 21 seconds
#########

## Parameters
expprepath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
windsize <- 200

## Read input tables
expdfpre <- read.csv(expprepath)

## Retrieving the bedgraph paths
bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
    package = "tepr"))

## Replace the path column of expdfpre with the previously retrieved paths
## and writing the new experiment file to the current folder
expdfpre$path <- bgpathvec
write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE, quote = FALSE)
expprepath <- "exptab-preprocessing.csv"

gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
genomename <- "hg38"

## The lines below are optional. The chromosome info can be retrieved automatically
## Make chromosome info retrieval explicit for building the vignette
chromtabtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
allchromvec <- GenomeInfoDb::seqnames(chromtabtest)
chromtabtest <- chromtabtest[allchromvec[which(allchromvec == "chr13")], ]

finaltab <- preprocessing(expprepath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename, finaltabpath = getwd(), finaltabname = "anno.tsv",
    chromtab = chromtabtest, showtime = FALSE, verbose = FALSE)


#########
# tepr analysis - takes ~ 1 seconds
#########

## Parameters (transpath limited to 6 transcripts)
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Read input tables
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)

reslist <- tepr(expdf, transdf, expthres, showtime = FALSE, verbose = FALSE)


## ----prequick, echo = FALSE---------------------------------------------------
print(head(finaltab, 3))


## ----teprquick, echo = FALSE--------------------------------------------------
message("The table meandifference:\n")
print(head(reslist[[1]], 2))

message("\n\n The table universegroup:\n")
print(head(reslist[[2]], 2))






## ----echo=FALSE, fig.cap="structure"------------------------------------------
knitr::include_graphics(system.file("extdata", "structure.png", package = "tepr"))


## ----retrieveanno-------------------------------------------------------------
allannobed <- retrieveanno(expprepath, gencodepath, verbose = FALSE)
message("\n The result is:\n")
print(head(allannobed, 3))


## ----makewindows--------------------------------------------------------------
allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)
message("\n The result is:\n")
print(head(allwindowsbed, 3))








## ----showfinaltab, echo = FALSE-----------------------------------------------
finaltabpath <- system.file("extdata", "finaltab-chr13.tsv", package = "tepr")
finaltab <- read.delim(finaltabpath, header = FALSE)
print(head(finaltab, 3))


## ----reading-anno-scores------------------------------------------------------
## Define manually the column names for display purpose
colnames(transdf) <- c("biotype", "chr", "coor1", "coor2", "transcript", "gene",
    "strand", "window", "id", "ctrl_rep1.plus", "ctrl_rep1.plus_score",
    "ctrl_rep1.minus", "ctrl_rep1.minus_score", "ctrl_rep2.plus",
    "ctrl_rep2.plus_score", "ctrl_rep2.minus", "ctrl_rep2.minus_score",
    "HS_rep1.plus", "HS_rep1.plus_score", "HS_rep1.minus",
    "HS_rep1.minus_score", "HS_rep2.plus", "HS_rep2.plus_score",
    "HS_rep2.minus", "HS_rep2.minus_score")

message("The table given by the preprocessing function is:\n")
print(head(transdf, 3))

message("\n The expdf table contains information about each replicate (here limited to one):\n")
head(expdf)


## ----checkexptab--------------------------------------------------------------
checkexptab(expdf)


## ----averageandfilterexprs----------------------------------------------------
resallexprs <- averageandfilterexprs(expdf, transdf, expthres, verbose = FALSE)


## ----exprsvec-----------------------------------------------------------------
print(resallexprs[[2]])


## ----countna------------------------------------------------------------------
rescountna <- countna(resallexprs, expdf, verbose = FALSE)
print(rescountna)


## ----ecdf---------------------------------------------------------------------
resecdflist <- genesECDF(resallexprs, expdf, verbose = FALSE)


## ----resecdf------------------------------------------------------------------
print(head(as.data.frame(resecdflist[[1]]), 3))


## ----nbwindows----------------------------------------------------------------
nbwindows <- resecdflist[[2]]
print(nbwindows)


## ----meandifference-----------------------------------------------------------
resmeandiff <- meandifference(resecdflist[[1]], expdf, nbwindows, verbose = FALSE)
print(head(resmeandiff, 3))


## ----splitmeans---------------------------------------------------------------
## Split the results by transcripts
bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))


## ----allauc-------------------------------------------------------------------
## Calculate Area Under Curve (AUC) and Differences of AUC for Transcript Data
resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
print(head(resauc, 3))


## ----knee---------------------------------------------------------------------
resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
print(resknee)


## ----attenuation--------------------------------------------------------------
resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
    resmeandiff, verbose = FALSE)
print(head(resatt, 3))














## ----universegroup------------------------------------------------------------
res <- universegroup(resatt, expdf, verbose = FALSE)
print(head(res, 2))


## ----plotecdf, warning = FALSE------------------------------------------------
colvec <- c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
plotecdf(resmeandiff, res, expdf, "EGFR", colvec, plot = TRUE, verbose = FALSE)


## ----plotaucgroups, warning = FALSE-------------------------------------------
plotauc(res, expdf, plottype = "groups", plot = TRUE)


## ----echo=FALSE, fig.cap="AUC group plot"-------------------------------------
knitr::include_graphics(system.file("extdata", "AUCcompare_group.png", package = "tepr"))


## ----plotaucpval, warning = FALSE---------------------------------------------
genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")
plotauc(res, expdf, genevec, plot = TRUE)


## ----echo=FALSE, fig.cap="AUC pval plot"--------------------------------------
knitr::include_graphics(system.file("extdata", "AUCcompare_pval.png", package = "tepr"))


## ----metageneplot, warning = FALSE--------------------------------------------
plotmetagenes(res, resmeandiff, expdf, plottype = "attenuation", plot = TRUE)


## ----echo=FALSE, fig.cap="Attenuation meta"-----------------------------------
knitr::include_graphics(system.file("extdata", "metagene_attenuation.png", package = "tepr"))


## ----kneepercent, warning = FALSE---------------------------------------------
plothistoknee(res, plot = TRUE)


## ----echo=FALSE, fig.cap="histo percent"--------------------------------------
knitr::include_graphics(system.file("extdata", "histo_percent.png", package = "tepr"))


## ----kneekb, warning = FALSE--------------------------------------------------
plothistoknee(res, plottype = "kb", plot = TRUE)


## ----echo=FALSE, fig.cap="histo kb"-------------------------------------------
knitr::include_graphics(system.file("extdata", "histo_kb.png", package = "tepr"))










## ----singlecond---------------------------------------------------------------
## The experiment table is limited to one condition and one replicate
expdfonecond <- expdf[which(expdf$condition == "HS" & expdf$replicate == 1), ]

## The table obtained by preprocessing is limited to the condition 'HS' and replicate 1
transdfonecond <- transdf[, c(seq_len(9), 18, 19, 20, 21)]

## Computing the object 'res' with tepr on one condition
resonecond <- tepr(expdfonecond, transdfonecond, expthres, controlcondname = "HS", verbose = FALSE)


## ----ecdfcond, warning = FALSE------------------------------------------------
plotecdf(resonecond[[1]], resonecond[[2]], expdfonecond, genename = "EGFR",
    colvec = c("#90AFBB"), plot = TRUE, verbose = FALSE)


## ----histocond, warning = FALSE-----------------------------------------------
## Randomly marking 3 transcripts as attenuated as a mock example
idxatt <- sample(seq_len(6), 3)
resonecond[[2]][idxatt, "Group"] <- "Attenuated"
resonecond[[2]][idxatt, "Universe"] <- TRUE
plothistoknee(resonecond[[2]], kneename = "knee_AUC_HS", plottype = "percent", plot = TRUE, verbose = FALSE)

