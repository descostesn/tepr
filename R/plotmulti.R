
## See also: teprmulti, plotecdf, plotauc, plotmetagenes

plotmulti <- function(resteprmulti, expdf, ecdfgenevec, genaucvec = NA,
    outfold = ".", digits = 2, middlewind = 100, pval = 0.01,
    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    aucaxisminx = -10, aucaxismaxx = 100, aucaxisminy = -10, aucaxismaxy = 100,
    aucmaintitle = "", aucsubtitle = "", auclegendpos = "bottom",
    formatname = "pdf", uniname = "Universe", groupname = "Group",
    histkneexlim = NA, binwidthvalhistknee = NA, verbose = TRUE) {

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

        ## Generate the plot of auc by groups
        if (verbose) message("\t ## plot auc by groups")
        expnamevec <- unique(expdf$condition)
        name1 <- expnamevec[1]
        name2 <- expnamevec[2]
        pvalks <- paste0("adjFDR_p_dAUC_Diff_meanFx_", name1, "_", name2)
        labelx <- paste0("AUC in ", name1)
        labely <- paste0("AUC in ", name2)
        aucfilename <- paste0("AUCcompare_groups_", name1, "_", name2)
        plotauc(tab = complist[[2]], genevec = genaucvec, auc_ctrlname = name1,
            auc_stressname = name2, pvalkstestcolname = pvalks, labelx = labelx,
            labely = labely, axismin_x = aucaxisminx, axismax_x = aucaxismaxx,
            axismin_y = aucaxisminy, axismax_y = aucaxismaxy,
            maintitle = aucmaintitle, subtitle = aucsubtitle,
            legendpos = auclegendpos, formatname = formatname,
            outfold = outfoldcomp, outfile = aucfilename, plottype = "groups",
            plot = FALSE, universename = uniname, groupname = groupname,
            verbose = verbose)

        ## Generate the plot of auc by pval
        if (!is.na(genaucvec)) {
            if (verbose) message("\t ## plot auc by pval for the given genes")
            aucfilename <- paste0("AUCcompare_pval_", name2, "_", name1)
            plotauc(tab = complist[[2]], genevec = genaucvec,
                auc_ctrlname = name1, auc_stressname = name2,
                pvalkstestcolname = pvalks, labelx = labelx, labely = labely,
                axismin_x = aucaxisminx, axismax_x = aucaxismaxx,
                axismin_y = aucaxisminy, axismax_y = aucaxismaxy,
                maintitle = aucmaintitle, subtitle = aucsubtitle,
                legendpos = auclegendpos, formatname = formatname,
                outfold = outfoldcomp, outfile = aucfilename, plottype = "pval",
                plot = FALSE, universename = uniname,
                groupname = groupname, verbose = verbose)
        }

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

        ## plothistoknee by percent
        if (verbose) message("\t ## plothistoknee by percent")
        kneename <- paste0("knee_AUC_", name2)
        plothistoknee(unigroupdf = complist[[2]], plottype = "percent",
            xlimvec = histkneexlim, binwidthval = binwidthvalhistknee, !!
            kneename = kneename, plot = FALSE, outfold = outfoldcomp,
            formatname = formatname = , universename = uniname,
            groupname = groupname, verbose = verbose)

        ## plothistoknee by kb

    }, resteprmulti, names(resteprmulti), MoreArgs = list(expdf, ecdfgenevec,
        genaucvec, colvec, digits, middlewind, pval, formatname, aucaxisminx,
        aucaxismaxx, aucaxisminy, aucaxismaxy, aucmaintitle, aucsubtitle,
        auclegendpos, uniname, groupname, histkneexlim,
        binwidthvalhistknee, verbose)))

}