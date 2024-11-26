retrieveanno <- function(saveobjectpath = NA) {

    if (!is.na(saveobjectpath) && !file.exists(saveobjectpath))
        dir.create(saveobjectpath, recursive = TRUE)
}