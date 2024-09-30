.callggplotauc <- function(df, aesvar, geompointinfo, geompointinfo2,
    geompointinfo3, plottype, axismin_x, axismax_x, axismin_y, axismax_y,
    labelx, labely, maintitle, subtitle, legendpos, plot, outfile, formatname,
    outfold, genevec) {

        ## Structure of the basic scatterplot
        g <- ggplot2::ggplot(df, aesvar) + geompointinfo + geompointinfo2

        if (isTRUE(all.equal(plottype, "pval")) && !is.na(genevec)) {

            ## Adding highlight of the genes
            g <- g + ggrepel::geom_label_repel(data = subset(df,
                gene %in% genevec), aes(label = gene), box.padding = 0.55, # nolint
                point.padding = 0, segment.color = "black", max.overlaps = 50,
                color = "red") +
                ggplot2::scale_color_gradient2(midpoint = 0, low = "white",
                    mid = "grey", high = "darkgreen")
        } else {
            g <- g + geompointinfo3
        }

        ## Formatting functions
        g <- g + ggplot2::xlim(axismin_x, axismax_x) +
            ggplot2::ylim(axismin_y, axismax_y) +
            ggplot2::labs(x = labelx, y = labely, legend = "-log10 p-value", # nolint
                color = "-log10 p-value", title = maintitle, # nolint
                subtitle = subtitle) +
            ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_classic() +
            ggplot2::theme(legend.position = legendpos)

        if (plot) {
            warning("You chose to plot the auc, the figure is not saved.") # nolint
            print(g)
        } else {
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}

.checkplotaucparams <- function(plottype, auc_ctrlname, auc_stressname,
    pvalkstestcolname, genevec, tab) {

        if (!isTRUE(all.equal(plottype, "pval")) &&
            !isTRUE(all.equal(plottype, "groups")))
                stop("plottype should be equal to 'pval' or 'groups'.")

        colnamevec <- c(auc_ctrlname, auc_stressname, pvalkstestcolname)
        .colnamecheck(colnamevec, tab)

        if (isTRUE(all.equal(plottype, "groups")) && !is.na(genevec[1]))
            stop("The vector of genes is not necessary for plotting groups")
}

plotauc <- function(tab, genevec = NA, # nolint
    auc_ctrlname = "AUC_ctrl", auc_stressname = "AUC_HS",
    pvalkstestcolname = "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in Stress", axismin_x = -10,
    axismax_x = 100, axismin_y = -10, axismax_y = 100, maintitle = "",
    subtitle = "", legendpos = "bottom", formatname = "pdf", outfold = "./",
    outfile = "AUCcompare_pval.pdf", plottype = "pval", plot = FALSE,
    universename = "Universe", groupname = "Group") {

        .checkplotaucparams(plottype, auc_ctrlname, auc_stressname,
            pvalkstestcolname, genevec, tab)

        if (isTRUE(all.equal(plottype, "pval"))) {
            df <- cbind(tab, kstestlog10 = -log10(tab[, pvalkstestcolname]))
            kstestlog10str <- "kstestlog10"
            df <- df %>% dplyr::arrange(df[, kstestlog10str])
            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname), # nolint
             color = !!sym(kstestlog10str))
            geompointinfo <- ggplot2::geom_point(size = 0.5)
            geompointinfo2 <- ggplot2::geom_density_2d()
        } else {
            df <- tab %>% dplyr::filter(!!sym(universename) == FALSE) # nolint
            dfatt <- tab %>% dplyr::filter(!!sym(groupname) == "Attenuated") # nolint
            dfoutgroup <- tab %>% dplyr::filter(!!sym(groupname) == "Outgroup") # nolint

            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname)) # nolint
            geompointinfo <- ggplot2::geom_point(size = 0.5, color = "grey")
            geompointinfo2 <- ggplot2::geom_point(data = dfatt, aesvar,
                color = "#e76f51", size = 1)
            geompointinfo3 <- ggplot2::geom_point(data = dfoutgroup, aesvar,
                    color = "#e9c46a", size = 1)
        }

        .callggplotauc(df, aesvar, geompointinfo, geompointinfo2,
            geompointinfo3, plottype, axismin_x, axismax_x, axismin_y,
            axismax_y, labelx, labely, maintitle, subtitle, legendpos, plot,
            outfile, formatname, outfold)
}
