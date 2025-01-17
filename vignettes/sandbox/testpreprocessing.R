library(tepr)

gencodepath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "./objects-drbttseq"
finaltabname <- "drbttseq.tsv"
saveobjectpath <- finaltabpath

windsize <- 200
blacklistshpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
nbcputrans <- 15
reload <- TRUE
verbose <- TRUE
showtime <- TRUE
showmemory <- TRUE
genomename <- "hg38"
tmpfold <- "./tmp-15cpu"
savefinaltable <- TRUE
deletetmp <- FALSE
chromtab <- NA

finaltable <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans, finaltabpath,
    finaltabname, tmpfold, saveobjectpath,
    savefinaltable, reload, showtime, showmemory,
    deletetmp, chromtab, verbose)
