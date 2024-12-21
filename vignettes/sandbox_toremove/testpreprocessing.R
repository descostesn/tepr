library("GenomeInfoDb")
library("GenomicRanges")
library("rtracklayer")
library("parallel")
library("purrr")
library("dplyr")
library("tidyverse")

gencodepath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
finaltabpath <- "./objects-tsv-7cpus"
finaltabname <- "cugusi.tsv"
saveobjectpath <- finaltabpath

windsize <- 200
blacklistshpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
nbcputrans <- 7
reload <- TRUE
verbose <- TRUE
showtime <- TRUE
showmemory <- TRUE
genomename <- "hg38"
tmpfold <- "./tmp-7cpu"
savefinaltable <- TRUE
deletetmp <- FALSE
# nbcpubg <- 1
# nbcputrans <- 1


#source("/g/romebioinfo/Projects/tepr/R/preprocessing-blacklisthighmap.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing-blacklisthighmap-v2.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing-blacklisthighmap-utils.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing-makewindows.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing-retrieveanno.R")
#source("/g/romebioinfo/Projects/tepr/R/preprocessing-createtablescores.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing-createtablescores-v2.R")
source("/g/romebioinfo/Projects/tepr/R/preprocessing.R")
source("/g/romebioinfo/Projects/tepr/R/utils.R")

finaltable <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans, finaltabpath,
    finaltabname, tmpfold, saveobjectpath,
    savefinaltable, reload, showtime, showmemory,
    deletetmp, verbose)
