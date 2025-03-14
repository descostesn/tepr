## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
    package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"


preprocessing <- function(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename)