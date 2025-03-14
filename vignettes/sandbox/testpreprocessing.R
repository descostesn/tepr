library("tepr")

gencodepath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/gencode.v43.basic.annotation.gtf" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
finaltabpath <- "objects-tsv-10cpus"
finaltabname <- "cugusi.tsv"
saveobjectpath <- finaltabpath

windsize <- 200
blacklistshpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/hg38-blacklist.v2.bed" # nolint
maptrackpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/k50.umap.hg38.0.8.bed" # nolint
nbcputrans <- 10
reload <- TRUE
verbose <- TRUE
showtime <- TRUE
showmemory <- TRUE
genomename <- "hg38"
tmpfold <- "tmp-10cpu"
savefinaltable <- TRUE
deletetmp <- FALSE
chromtest <- rtracklayer::SeqinfoForUCSCGenome(genomename)
idxkeep <- GenomeInfoDb::seqnames(chromtest)[grep("_|chrM",
        GenomeInfoDb::seqnames(chromtest), perl = TRUE, invert = TRUE)]
chromtest <- chromtest[idxkeep,]


finaltable <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistshpath, genomename, nbcputrans, finaltabpath,
    finaltabname, tmpfold, saveobjectpath,
    savefinaltable, reload, showtime, showmemory,
    deletetmp, chromtab = chromtest, verbose)

