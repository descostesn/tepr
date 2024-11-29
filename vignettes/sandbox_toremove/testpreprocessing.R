library("GenomeInfoDb")
library("GenomicRanges")
library("rtracklayer")
library("parallel")
library("purrr")
library("dplyr")


gencodepath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing"
finaltabname <- "cugusi.tsv"
saveobjectpath <- finaltabpath

windsize <- 200
blacklistshpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
nbcpubg <- 8
nbcputrans <- 20

finaltab <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, nbcputrans, nbcpubg, finaltabpath, finaltabname,
    saveobjectpath)
