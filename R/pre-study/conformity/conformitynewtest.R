library("rtracklayer")


bgvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/protein_coding_score/ctrl_rep1.forward.window200.MANE.wmean.name.score"

allbgnicpath <- "/g/romebioinfo/tmp/preprocessing/backup/bedgraphwmeanlist.rds"

bgvic <- read.delim(bgvicpath, header = FALSE)
allbgnic <- readRDS(allbgnicpath)
