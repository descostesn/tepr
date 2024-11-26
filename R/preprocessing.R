retrieveanno <- function(exptabpath, gencodepath, saveobjectpath = NA,
    verbose = TRUE) {

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Reading the information about experiments
    if (verbose) message("Reading the information about experiments")
    exptab <- read.csv(exptabpath, header = TRUE)
    checkexptab(exptab) # nolint

    ## Reading gencode file
    if (verbose) message("Reading gencode file")
    gencode <- read.delim(gencodepath, header = FALSE, skip = 5)

}