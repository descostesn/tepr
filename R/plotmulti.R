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

.multiplotauc <- function(nameone, nametwo, complist, genaucvec, aucaxisminx,
    aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
    auclegendpos, formatname, outfoldcomp, uniname, groupname, verbose) {

        ## Generate the plot of auc by groups
        if (verbose) message("\t ## plot auc by groups")
        pvalks <- paste0("adjFDR_p_dAUC_Diff_meanFx_", nametwo, "_",
            nameone)
        labelx <- paste0("AUC in ", nameone)
        labely <- paste0("AUC in ", nametwo)
        aucfilename <- paste0("AUCcompare_groups_", nameone, "_",
            nametwo)
        plotauc(tab = complist[[2]],
            auc_ctrlname = paste0("AUC_", nameone),
            auc_stressname = paste0("AUC_", nametwo),
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
        if (!is.na(genaucvec[1])) {
            if (verbose) message("\t ## plot auc by pval for the ",
                "given genes")
            aucfilename <- paste0("AUCcompare_pval_", nametwo, "_", nameone)
            plotauc(tab = complist[[2]], genevec = genaucvec,
                auc_ctrlname = paste0("AUC_", nameone),
                auc_stressname = paste0("AUC_", nametwo),
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

.multiplotmetagenes <- function(complist, nameone, nametwo, formatname,
    outfoldcomp, verbose) {

        daucname <- paste0("dAUC_Diff_meanFx_", nametwo, "_", nameone)
        aucctrlname <- paste0("AUC_", nameone)
        aucstressname <- paste0("AUC_", nametwo)

        ## Plot metagene by attenuation
        if (verbose) message("\t ## Plot metagene by attenuation")
        plotmetagenes(unigroupdf = complist[[2]], dfmeandiff = complist[[1]],
            plottype = "attenuation", daucname = daucname,
            auc_ctrlname = aucctrlname, auc_stressname = aucstressname,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by outgroup
        if (verbose) message("\t ## Plot metagene by outgroup")
        plotmetagenes(unigroupdf = complist[[2]], dfmeandiff = complist[[1]],
            plottype = "outgroup", daucname = daucname,
            auc_ctrlname = aucctrlname, auc_stressname = aucstressname,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by universe
        if (verbose) message("\t ## Plot metagene by universe")
        plotmetagenes(unigroupdf = complist[[2]], dfmeandiff = complist[[1]],
            plottype = "universe", daucname = daucname,
            auc_ctrlname = aucctrlname, auc_stressname = aucstressname,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)

        ## Plot metagene by all
        if (verbose) message("\t ## Plot metagene for all transcripts")
        plotmetagenes(unigroupdf = complist[[2]], dfmeandiff = complist[[1]],
            plottype = "all", daucname = daucname,
            auc_ctrlname = aucctrlname, auc_stressname = aucstressname,
            plot = FALSE, formatname = formatname, outfold = outfoldcomp,
            verbose = verbose)
}

.multiplothistoknee <- function(complist, histkneexlim, binwidthvalhistknee,
            outfoldcomp, formatname, uniname, groupname, nametwo, verbose) {

                kneename <- paste0("knee_AUC_", nametwo)

                ## plothistoknee by percent
                if (verbose) message("\t ## plothistoknee by percent")
                plothistoknee(unigroupdf = complist[[2]], plottype = "percent",
                    xlimvec = histkneexlim, binwidthval = binwidthvalhistknee,
                    kneename = kneename, plot = FALSE, outfold = outfoldcomp,
                    formatname = formatname, universename = uniname,
                    groupname = groupname, verbose = verbose)

                ## plothistoknee by kb
                if (verbose) message("\t ## plothistoknee by kb")
                plothistoknee(unigroupdf = complist[[2]], plottype = "kb",
                    xlimvec = histkneexlim, binwidthval = binwidthvalhistknee,
                    kneename = kneename, plot = FALSE, outfold = outfoldcomp,
                    formatname = formatname, universename = uniname,
                    groupname = groupname, verbose = verbose)
}

#' Generate all tepr plots for all experiment comparisons
#'
#' @description
#' This function generates for all experiment comparisons contained in the
#' object \code{resteprmulti} all plots of tepr: ECDF, auc, metagene, and
#' histtoknee.
#'
#' @usage
#' plotmulti(resteprmulti, expdf, ecdfgenevec, outfold = tempdir(), digits = 2,
#' middlewind = 100, pval = 0.01, colvec = c("#90AFBB", "#10AFBB",
#' "#FF9A04", "#FC4E07"), genaucvec = NA, aucaxisminx = -10,
#' aucaxismaxx = 100, aucaxisminy = -10, aucaxismaxy = 100, aucmaintitle = "",
#' aucsubtitle = "", auclegendpos = "bottom", formatname = "pdf",
#' uniname = "Universe", groupname = "Group", histkneexlim = NA,
#' binwidthvalhistknee = NA, verbose = TRUE)
#' 
#' @param resteprmulti Result returned by the function \code{teprmulti}.
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param ecdfgenevec A vector specifying the names of the genes of interest to
#'  plot the ecdf of.
#' @param outfold Path to the output folder where the plots will be written.
#'  Subfolders with the names of the comparisons are automatically created.
#'  Default is \code{tempdir()}.
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
#' @param formatname Format of the saved plot (e.g., "pdf", "png"). Default is
#'  \code{"pdf"}.
#' @param uniname Column name in the second element of \code{resteprmulti}
#'  representing the universe selection. Default is \code{"Universe"}.
#' @param groupname Column name in the second element of \code{resteprmulti}
#'  representing the type of group a transcript belong to. Default is
#' \code{"Group"}.
#' @param histkneexlim For the plot histoknee, a numeric vector of length 2
#'  specifying the limits of the x-axis. Default is \code{NA}, which
#' automatically sets the limits based on \code{plottype}.
#' @param binwidthvalhistknee For the plot histoknee, a numeric value for the
#'  width of the bins in the histogram. Default is \code{NA}, which
#'  automatically selects a bin width based on \code{plottype}.
#' @param verbose A logical flag indicating whether to display detailed
#'  messages about the function's progress. Default is \code{TRUE}.
#'
#' @return Nothing is returned. Figures are written to outfold in the subfolder
#'  of the corresponding comparison.
#'
#' @details
#' The function goes through each element of resteprmulti which corresponds to
#' a comparison of two conditions. For each element it calls the following
#' functions:
#' #' \itemize{
#'   \item \code{"plotecdf"}: The function generates a figure for each gene
#'  given in ecdfgenevec.
#'   \item \code{"plotauc"}: Generates figures by groups and pval. The lattest
#'  figure is not generated if genaucvec = NA.
#'   \item \code{"plotmetagenes"}: Generates the figures by \code{attenuation},
#'  \code{outgroup}, \code{universe}, and \code{all}.
#'   \item \code{"plothistoknee"}: Generate the figures by \code{percent} and
#'  \code{kb}.
#' }
#'
#' @examples
#' # Assuming resteprmulti is the object returned by the function teprmulti
#' # and expdf contains the necessary data:
#' \dontrun{
#' plotmulti(resteprmulti, expdf, ecdfgenevec = c("EGFR", "DAP", "FLI1"))}
#'
#' @seealso
#' [teprmulti], [plotecdf], [plotauc], [plotmetagenes], [plothistoknee]
#'
#' @export

plotmulti <- function(resteprmulti, expdf, ecdfgenevec, outfold = tempdir(),
    digits = 2, middlewind = 100, pval = 0.01,
    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    genaucvec = NA, aucaxisminx = -10, aucaxismaxx = 100, aucaxisminy = -10,
    aucaxismaxy = 100, aucmaintitle = "", aucsubtitle = "",
    auclegendpos = "bottom", formatname = "pdf", uniname = "Universe",
    groupname = "Group", histkneexlim = NA, binwidthvalhistknee = NA,
    verbose = TRUE) {

    if (!length(unique(expdf$condition)) > 2)
        stop("\n[tepr] Error: Too few conditions.\n",
            "  plotmulti requires teprmulti output (>2 conditions).\n")

    invisible(mapply(function(complist, compname, expdf, ecdfgenevec,
        genaucvec, colvec, digits, middlewind, pval, formatname, aucaxisminx,
        aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
        auclegendpos, uniname, groupname, histkneexlim, binwidthvalhistknee,
        verbose) {

        if (verbose) message("\n Generating plots for ", compname)
        outfoldcomp <- file.path(outfold, compname)
        expnamevec <- unlist(strsplit(compname, "_vs_"))
        nameone <- expnamevec[1]
        nametwo <- expnamevec[2]
        idxtwocond <- which(expdf$condition == nameone | expdf$condition == nametwo)
        expdftwocond <- expdf[idxtwocond, ]

        ## Generating the plot of the ecdf empirical distribution and
        ## nsc-rna-seq signal
        .multiplotecdf(ecdfgenevec, complist, expdftwocond, colvec, outfoldcomp,
            digits, middlewind, pval, formatname, verbose)

        ## Generate the plot of auc by groups and pval
        .multiplotauc(nameone, nametwo, complist, genaucvec, aucaxisminx,
            aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
            auclegendpos, formatname, outfoldcomp, uniname, groupname, verbose)

        ## Plot metagene by attenuation, outgroup, universe, and all
        .multiplotmetagenes(complist, nameone, nametwo, formatname, outfoldcomp,
            verbose)

        ## plothistoknee by percent and kb
        .multiplothistoknee(complist, histkneexlim, binwidthvalhistknee,
            outfoldcomp, formatname, uniname, groupname, nametwo, verbose)

    }, resteprmulti, names(resteprmulti), MoreArgs = list(expdf, ecdfgenevec,
        genaucvec, colvec, digits, middlewind, pval, formatname, aucaxisminx,
        aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
        auclegendpos, uniname, groupname, histkneexlim,
        binwidthvalhistknee, verbose)))

}
