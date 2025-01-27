
## See also: teprmulti, plotecdf

plotmulti <- function(resteprmulti, expdf, ecdfgenevec, outfold = ".",
    digits = 2, colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    middlewind = 100, pval = 0.01, formatname = "pdf", verbose = TRUE) {

    if (!length(unique(expdf$condition)) > 2)
        stop("There are less than two conditions in your experiment ",
            "table. The input list must be the result of teprmulti.")

    ## complist <- resteprmulti[[1]]; compname <- names(resteprmulti)[1]
    invisible(mapply(function(complist, compname, expdf, ecdfgenevec, colvec,
        digits, middlewind, pval, formatname, verbose) {

        if (verbose) message("Generating plots for ", compname)
        outfoldcomp <- file.path(outfold, compname)

        ## Generating the plot of the ecdf empirical distribution and
        ## nsc-rna-seq signal
        if (verbose) message("\t ## plot ecdf")
        sapply(ecdfgenevec, function(currentgene, complist, expdf, colvec,
            outfoldcomp, digits, middlewind, pval, formatname, verbose) {
            if (verbose) message("\t\t plot ecdf for ", currentgene)
            plotecdf(dfmeandiff = complist[[1]], unigroupdf = complist[[2]],
                expdf = expdf, genename = currentgene, colvec = colvec,
                outfold = outfoldcomp, digits = digits, middlewind = middlewind,
                pval = pval, plot = FALSE, formatname = formatname,
                verbose = verbose)
        }, complist, expdf, colvec, outfoldcomp, digits, middlewind,
                pval, formatname, verbose)

    }, resteprmulti, names(resteprmulti), MoreArgs = list(expdf, ecdfgenevec,
        colvec, digits, middlewind, pval, formatname, verbose)))

}