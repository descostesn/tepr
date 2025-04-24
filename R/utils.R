.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("\n\t The element ", xname, # nolint
                " should be unique, contact the developer.\n") # nolint
}

.extractstr <- function(transtable) {

    str <- as.character(unique(transtable$strand))
    .checkunique(str, "str")
    if (isTRUE(all.equal(str, "+"))) {
        str <- "plus"
    } else if (isTRUE(all.equal(str, "-"))) {
        str <- "minus"
    } else {
        stop("\n\t In .computeecdf or countna, strand is neither plus or",
            " minus in the table returned by the function ",
            "averageandfilterexprs. This should not happen. Contact the ",
            "developer.\n")
    }
    return(str)
}

.colnamecheck <- function(colnamevec, tab) {
            invisible(sapply(colnamevec, function(currentcol, tab) {
            idx <- grep(currentcol, colnames(tab))
            if (isTRUE(all.equal(length(idx), 0)))
                stop("\n\t The column ", currentcol, " does not exist in the ",
                    "provided table.\n")
        }, tab))
}

.buildcolnames <- function(expdf, alldf) {

    infocolnames <- c("biotype", "chr", "coor1", "coor2", "transcript",
        "gene", "strand", "window", "id")
    expcolnames <- unlist(apply(expdf, 1, function(x) {
        res <- paste0(x["condition"], "_rep", x["replicate"], ".", x["strand"])
        return(c(res, paste(res, "score", sep = "_")))
    }, simplify = FALSE))
    colnames(alldf) <- c(infocolnames, expcolnames)
    return(alldf)
}

.returnexpcolnames <- function(expdf) {
    expcolnames <- unlist(apply(expdf, 1, function(x) {
        res <- paste0(x["condition"], "_rep", x["replicate"], ".", x["strand"])
        return(c(res, paste(res, "score", sep = "_")))
    }, simplify = FALSE))
    return(expcolnames)
}

.dontcompare <- function(dontcompare, expdf, verbose) {

    ## Retrieve the condition names without duplicates
    condvec <- unique(expdf$condition)

    ## Create matrix with all comparisons
    matcond <- utils::combn(condvec, 2, simplify = TRUE)

    if (!is.null(dontcompare)) {

        if (!is.vector(dontcompare))
            stop("\n The variable dontcompare should be a vector.\n")

        ## Building all comparisons from matcond
        compvec <- apply(matcond, 2, function(x) paste0(x[1], "_vs_", x[2]))

        ## Retrieving the comparisons to exclude
        idx <- match(dontcompare, compvec)
        idxna <- which(is.na(idx))

        if (isTRUE(all.equal(length(idx), 0)) ||
            !isTRUE(all.equal(length(idxna), 0)))
            stop("\n Problem with the values contained in the dontcompare ",
                "vector. Make sure that your vector contains one of these:\n",
                paste(compvec, collapse = " - "))

        if (isTRUE(all.equal(length(dontcompare), ncol(matcond))))
            stop("\n All comparisons are removed, the function cannot be ",
                "executed\n")

        matcond <- matcond[, -idx]

        if (verbose) message("The following comparisons were excluded:\n ",
            paste(dontcompare, collapse = " - "))
    }
    return(matcond)
}



#' Check Validity of Experiment Table
#'
#' @description
#' The `checkexptab` function verifies the structure and content of an
#' experiment table to ensure it meets specific formatting requirements. It
#' checks for the presence of required columns, and validates that the
#' `direction` and `strand` columns contain only allowable values.
#'
#' @usage
#' checkexptab(exptab)
#'
#' @param exptab A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#'
#' @return
#' If the experiment table is valid, the function returns `NULL`. If the table
#' is invalid, the function throws an error specifying the issue.
#'
#' @details
#' The function performs the following checks:
#' - The column names of `exptab` must match exactly: `"condition"`,
#'  `"replicate"`, `"direction"`, and `"strand"`.
#' - The `direction` column must contain only `"forward"` and `"reverse"`.
#' - The `strand` column must contain only `"plus"` and `"minus"`.
#'
#' @examples
#' # Create a valid experiment table
#' exptab <- data.frame(
#'   condition = c("cond1", "cond2"),
#'   replicate = c(1, 1),
#'   direction = c("forward", "reverse"),
#'   strand = c("plus", "minus"),
#'   path = c("toto/", "toto/"))
#' checkexptab(exptab)  # Should pass without errors
#'
#' # Invalid experiment table (wrong column names)
#' invalid_exptab <- data.frame(
#'     cond = c("cond1", "cond2"),
#'     rep = c(1, 1),
#'     dir = c("forward", "reverse"),
#'     str = c("+", "-"),
#'     paths = c("toto/", "toto/"))
#' try(checkexptab(invalid_exptab))
#'
#' @export

checkexptab <- function(exptab) {

    colnamevec <- c("condition", "replicate", "direction", "strand", "path")
    if (!isTRUE(all.equal(sort(colnames(exptab)), sort(colnamevec))))
        stop("\n\t The experiment table should have the columns: ",
            "'condition', 'replicate', 'direction', 'strand', 'path'.\n")

    directionvec <- unique(exptab$direction)
    if (!isTRUE(all.equal(length(directionvec), 2)) ||
        !isTRUE(all.equal(directionvec, c("forward", "reverse"))))
        stop("\n\t Only two values are allowed for the column direction of the",
            "experiment table, 'forward' and 'reverse'.\n")

    strandvec <- unique(exptab$strand)
    if (!isTRUE(all.equal(strandvec, c("plus", "minus"))))
        stop("\n\t The strand column of the experiment table should only ",
            "contain 'plus' and 'minus'.\n")

    idxchar <- grep("_|-", exptab$condition)
    if (!isTRUE(all.equal(length(idxchar), 0)))
        stop("\n\t The condition names should not contain any special ",
            "characters such as '_' or '-'.\n")
}



#' Retrieve all the comparison names from the experiment table
#'
#' @description
#' The `showallcomp` function build the string of all comparisons possible
#' using the condition column of a provided experiment table.
#'
#' @usage
#' showallcomp(expdf, verbose = FALSE)
#'
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param verbose A logical flag indicating whether to print progress messages.
#'                Defaults to \code{FALSE}.
#'
#' @return
#' If less than three conditions, nothing. Otherwise a character vector of all
#' comparisons.
#'
#' @examples
#' # Create a valid experiment table
#' exptab <- data.frame(
#'   condition = c("cond1", "cond2", "cond3"),
#'   replicate = c(1, 1, 1),
#'   direction = c("forward", "reverse", "forward"),
#'   strand = c("plus", "minus", "plus"),
#'   path = c("toto/", "toto/", "toto/"))
#' checkexptab(exptab)
#' showallcomp(exptab)
#'
#' @importFrom utils combn
#' @export

showallcomp <- function(expdf, verbose = FALSE) {

    condvec <- unique(expdf$condition)

    if (isTRUE(all.equal(length(condvec), 2))) {
        message("\n Your table has only two conditions: ", condvec[1], "_vs_",
            condvec[2])
    } else if (isTRUE(all.equal(length(condvec), 1))) {
        message("\n Your table has only one conditions: ", condvec)
    } else {

        matcond <- utils::combn(condvec, 2, simplify = TRUE)
        compvec <- apply(matcond, 2, function(x) paste0(x[1], "_vs_", x[2]))
        if (verbose) message("All comparisons: ", paste(compvec,
            collapse = " - "))
        return(compvec)
    }
}
