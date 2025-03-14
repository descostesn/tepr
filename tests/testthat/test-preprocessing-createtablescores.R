## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"
chromtabtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
allchromvec <- GenomeInfoDb::seqnames(chromtabtest)
chromtabtest <- chromtabtest[allchromvec[which(allchromvec == "chr13")], ]

## Necessary result to call createtablescores
allannobed <- retrieveanno(exptabpath, gencodepath, verbose = FALSE)
allwindowsbed <- makewindows(allannobed, windsize, verbose = FALSE)
blacklisthighmap(maptrackpath, blacklistpath, exptabpath, nbcputrans = 1,
    allwindowsbed, windsize, genomename = genomename, chromtab = chromtabtest,
    verbose = FALSE)

## Calling the function to test
finaltabtest <- createtablescores(tmpfold = file.path(getwd(), "tmptepr"),
    exptabpath, savefinaltable = FALSE, verbose = FALSE)
