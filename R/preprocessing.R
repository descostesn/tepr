retrieveanno <- function(exptabpath, saveobjectpath = NA) {

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)

    ## Reading the information about experiments
    exptab <- read.csv(exptabpath, header = TRUE)
    checkexptab(exptab)

}