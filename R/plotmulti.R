
## See also: teprmulti

plotmulti <- function(resteprmulti, expdf, genename, outfold = ".", digits = 2,
    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    middlewind = 100, pval = 0.01, formatname = "pdf", verbose = TRUE) {

    if (!length(unique(expdf$condition)) > 2)
        stop("There are less than two conditions in your experiment ",
            "table. The input list must be the result of teprmulti.")

    ## complist <- resteprmulti[[1]]; compname <- names(resteprmulti)[1]
    invisible(mapply(function(complist, compname, expdf, genename, colvec,
        digits, middlewind, pval, formatname, verbose) {

        if (verbose) message("Generating plots for ", compname)
        outfoldcomp <- file.path(outfold, compname)

        ## Generating the plot of the ecdf empirical distribution and
        ## nsc-rna-seq signal
        plotecdf(dfmeandiff = complist[[1]], unigroupdf = complist[[2]],
            expdf = expdf, genename = genename, colvec = colvec,
            outfold = outfoldcomp, digits = digits, middlewind = middlewind,
            pval = pval, plot = FALSE, formatname = formatname,
            verbose = verbose)

    }, resteprmulti, names(resteprmulti), MoreArgs = list(expdf, genename,
        colvec, digits, middlewind, pval, formatname, verbose)))

}