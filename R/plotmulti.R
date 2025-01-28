#' Generate all tepr plots for all experiment comparisons
#'
#' @description
#' This function generates for all experiment comparisons contained in the
#' object \code{resteprmulti} all plots of tepr: ECDF, auc, metagene, and
#' histtoknee.
#'
#' @usage
#' plotmulti(resteprmulti, expdf, ecdfgenevec, outfold = ".", digits = 2,
#' middlewind = 100, pval = 0.01, colvec = c("#90AFBB", "#10AFBB", "#FF9A04",
#' "#FC4E07"), aucaxisminx = -10, aucaxismaxx = 100, aucaxisminy = -10,
#' aucaxismaxy = 100, aucmaintitle = "", aucsubtitle = "",
#' auclegendpos = "bottom", formatname = "pdf", uniname = "Universe",
#' groupname = "Group", histkneexlim = NA, binwidthvalhistknee = NA,
#' verbose = TRUE)
#'
#' @param resteprmulti Result returned by the function \code{teprmulti}.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param ecdfgenevec A vector specifying the names of the genes of interest to
#'  plot the ecdf of.
#' @param outfold Path to the output folder where the plots will be written.
#'  Subfolders with the names of the comparisons are automatically created.
#' @param digits For the ecdf plot, the number of decimal places to round the
#'  AUC and KS values. Default is \code{2}.
#' @param middlewind For the ecdf plot, the index of the middle window
#'  representing the region centered around the TSS. Default is \code{100}.
#' @param pval For the ecdf plot, a numeric value for the p-value threshold to
#'  determine the significance of the KS test. Default is \code{0.01}.
#' @param colvec For the ecdf plot, a vector of 4 colors used to distinguish
#'  the different conditions. Default is \code{c("#90AFBB", "#10AFBB",
#'  "#FF9A04", "#FC4E07")}.
#' @param genaucvec For the auc plot, vector of gene names to highlight,
#'  Used for the plot of type "pval". Default is \code{NA}. If left to NA, the
#'  plot taking into account the p-values is not generated.
#' @param aucaxisminx For the auc plot, minimum value for the x-axis. Default
#'  is \code{-10}.
#' @param aucaxismaxx For the auc plot, maximum value for the x-axis. Default
#'  is \code{100}.
#' @param aucaxisminy For the auc plot, minimum value for the y-axis. Default
#'  is \code{-10}.
#' @param aucaxismaxy For the auc plot, maximum value for the y-axis. Default
#'  is \code{100}.
#' @param aucmaintitle For the auc plot, main title of the plot. Default is an
#' empty string.
#' @param aucsubtitle For the auc plot, subtitle of the plot. Default is an
#' empty string.
#' @param auclegendpos For the auc plot, position of the legend. Default is
#' \code{"bottom"}.
#' 
#' 
#' formatname = "pdf", uniname = "Universe",
#' groupname = "Group", histkneexlim = NA, binwidthvalhistknee = NA,
#' verbose = TRUE)








## See also: teprmulti, plotecdf, plotauc, plotmetagenes, plothistoknee

.multiplotecdf <- function(ecdfgenevec, complist, expdf, colvec, outfoldcomp,
    digits, middlewind, pval, formatname, verbose) {

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
}

.multiplotauc <- function(name1, name2, complist, genaucvec, pvalks, labelx,
            labely, aucaxisminx, aucaxismaxx, aucaxisminy, aucaxismaxy,
            aucmaintitle, aucsubtitle, auclegendpos, formatname,
            outfoldcomp, aucfilename, uniname, groupname, verbose) {

                ## Generate the plot of auc by groups
                if (verbose) message("\t ## plot auc by groups")
                pvalks <- paste0("adjFDR_p_dAUC_Diff_meanFx_", name1, "_",
                    name2)
                labelx <- paste0("AUC in ", name1)
                labely <- paste0("AUC in ", name2)
                aucfilename <- paste0("AUCcompare_groups_", name1, "_",
                    name2)
                plotauc(tab = complist[[2]], genevec = genaucvec,
                    auc_ctrlname = name1, auc_stressname = name2,
                    pvalkstestcolname = pvalks, labelx = labelx,
                    labely = labely, axismin_x = aucaxisminx,
                    axismax_x = aucaxismaxx, axismin_y = aucaxisminy,
                    axismax_y = aucaxismaxy, maintitle = aucmaintitle,
                    subtitle = aucsubtitle, legendpos = auclegendpos,
                    formatname = formatname, outfold = outfoldcomp,
                    outfile = aucfilename, plottype = "groups",
                    plot = FALSE, universename = uniname, groupname = groupname,
                    verbose = verbose)

                ## Generate the plot of auc by pval
                if (!is.na(genaucvec)) {
                    if (verbose) message("\t ## plot auc by pval for the ",
                        "given genes")
                    aucfilename <- paste0("AUCcompare_pval_", name2, "_", name1)
                    plotauc(tab = complist[[2]], genevec = genaucvec,
                        auc_ctrlname = name1, auc_stressname = name2,
                        pvalkstestcolname = pvalks, labelx = labelx,
                        labely = labely, axismin_x = aucaxisminx,
                        axismax_x = aucaxismaxx, axismin_y = aucaxisminy,
                        axismax_y = aucaxismaxy, maintitle = aucmaintitle,
                        subtitle = aucsubtitle, legendpos = auclegendpos,
                        formatname = formatname, outfold = outfoldcomp,
                        outfile = aucfilename, plottype = "pval",
                        plot = FALSE, universename = uniname,
                        groupname = groupname, verbose = verbose)
                }
}

.multiplotmetagenes <- function(complist, daucname, name1, name2, formatname,
    outfoldcomp, verbose) {

        ## Plot metagene by attenuation
        if (verbose) message("\t ## Plot metagene by attenuation")
        daucname <- paste0("dAUC_Diff_meanFx_", name2, "_", name1)
        plotmetagenes(unigroupdf = complist[[2]],
            dfmeandiff = complist[[2]], plottype = "attenuation",
            daucname = daucname, auc_ctrlname = name1, auc_stressname = name2,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by outgroup
        if (verbose) message("\t ## Plot metagene by outgroup")
        plotmetagenes(unigroupdf = complist[[2]],
            dfmeandiff = complist[[2]], plottype = "outgroup",
            daucname = daucname, auc_ctrlname = name1, auc_stressname = name2,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by universe
        if (verbose) message("\t ## Plot metagene by universe")
        plotmetagenes(unigroupdf = complist[[2]],
            dfmeandiff = complist[[2]], plottype = "universe",
            daucname = daucname, auc_ctrlname = name1, auc_stressname = name2,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by all
        if (verbose) message("\t ## Plot metagene for all transcripts")
        plotmetagenes(unigroupdf = complist[[2]],
            dfmeandiff = complist[[2]], plottype = "all",
            daucname = daucname, auc_ctrlname = name1, auc_stressname = name2,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)
}

.multiplothistoknee <- function(complist, histkneexlim, binwidthvalhistknee,
            kneename, outfoldcomp, formatname, uniname, groupname, name2,
            verbose) {

                ## plothistoknee by percent
                if (verbose) message("\t ## plothistoknee by percent")
                kneename <- paste0("knee_AUC_", name2)
                plothistoknee(unigroupdf = complist[[2]], plottype = "percent",
                    xlimvec = histkneexlim, binwidthval = binwidthvalhistknee,
                    kneename = kneename, plot = FALSE, outfold = outfoldcomp,
                    formatname = formatname, universename = uniname,
                    groupname = groupname, verbose = verbose)

                ## plothistoknee by kb
                if (verbose) message("\t ## plothistoknee by kb")
                kneename <- paste0("knee_AUC_", name2)
                plothistoknee(unigroupdf = complist[[2]], plottype = "kb",
                    xlimvec = histkneexlim, binwidthval = binwidthvalhistknee,
                    kneename = kneename, plot = FALSE, outfold = outfoldcomp,
                    formatname = formatname, universename = uniname,
                    groupname = groupname, verbose = verbose)
}

plotmulti <- function(resteprmulti, expdf, ecdfgenevec, outfold = ".",
    digits = 2, middlewind = 100, pval = 0.01,
    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    genaucvec = NA, aucaxisminx = -10, aucaxismaxx = 100, aucaxisminy = -10,
    aucaxismaxy = 100, aucmaintitle = "", aucsubtitle = "",
    auclegendpos = "bottom", formatname = "pdf", uniname = "Universe",
    groupname = "Group", histkneexlim = NA, binwidthvalhistknee = NA,
    verbose = TRUE) {

    if (!length(unique(expdf$condition)) > 2)
        stop("There are less than two conditions in your experiment ",
            "table. The input list must be the result of teprmulti.")

    ## complist <- resteprmulti[[1]]; compname <- names(resteprmulti)[1]
    invisible(mapply(function(complist, compname, expdf, ecdfgenevec,
        genaucvec, colvec, digits, middlewind, pval, formatname, aucaxisminx,
        aucaxismaxx, aucaxisminy, aucaxismaxy, aucsubtitle, auclegendpos,
        uniname, groupname, histkneexlim, binwidthvalhistknee, verbose) {

        if (verbose) message("\n Generating plots for ", compname)
        outfoldcomp <- file.path(outfold, compname)
        expnamevec <- unique(expdf$condition)
        name1 <- expnamevec[1]
        name2 <- expnamevec[2]

        ## Generating the plot of the ecdf empirical distribution and
        ## nsc-rna-seq signal
        .multiplotecdf(ecdfgenevec, complist, expdf, colvec, outfoldcomp,
            digits, middlewind, pval, formatname, verbose)

        ## Generate the plot of auc by groups and pval
        .multiplotauc(name1, name2, complist, genaucvec, pvalks, labelx, labely,
            aucaxisminx, aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle,
            aucsubtitle, auclegendpos, formatname, outfoldcomp, aucfilename,
            uniname, groupname, verbose)

        ## Plot metagene by attenuation, outgroup, universe, and all
        .multiplotmetagenes(complist, daucname, name1, name2, formatname,
            outfoldcomp, verbose)

        ## plothistoknee by percent and kb
        .multiplothistoknee(complist, histkneexlim, binwidthvalhistknee,
            kneename, outfoldcomp, formatname, uniname, groupname, verbose)

    }, resteprmulti, names(resteprmulti), MoreArgs = list(expdf, ecdfgenevec,
        genaucvec, colvec, digits, middlewind, pval, formatname, aucaxisminx,
        aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
        auclegendpos, uniname, groupname, histkneexlim,
        binwidthvalhistknee, verbose)))

}
