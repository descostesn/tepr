plothistoknee <- function(unigroupdf, plottype = "percent", xlimvec = NA, # nolint
    binwidthval = NA, kneename = "knee_AUC_HS", plot = FALSE, outfold = "./",
    formatname = "pdf", universename = "Universe", groupname = "Group") {

        if (!isTRUE(all.equal(plottype, "percent")) &&
            !isTRUE(all.equal(plottype, "kb")))
                stop("Plot type should be percent or kb")

        colnamevec <- c(universename, groupname, kneename)
        .colnamecheck(colnamevec, unigroupdf)

        if (isTRUE(all.equal(plottype, "percent"))) {
            gtypeaes <- ggplot2::aes(x = !!sym(kneename) / 2)
            if (is.na(xlimvec)) xlimvec <- c(0, 100)
            if (is.na(binwidthval)) binwidthval <- 5
            gtypelabs <- ggplot2::xlab("Distance TSS to knee (% of the gene)")
        } else {
            gtypeaes <- ggplot2::aes(x = (!!sym(kneename) * window_size) / 1000)
            if (is.na(xlimvec)) xlimvec <- c(0, 350)
            if (is.na(binwidthval)) binwidthval <- 10
            gtypelabs <- ggplot2::labs(x = "Distance TSS to knee (kb)",
                y = "Count",
                title = "TSS to knee position in kb for attenuated genes")
        }

        g <- ggplot2::ggplot(unigroupdf %>%
            dplyr::filter(!!sym(universename) == TRUE &
                !!sym(groupname) == "Attenuated"), gtypeaes) +
            ggplot2::geom_histogram(binwidth = binwidthval, fill = "grey",
                color = "black", boundary = 0) + ggplot2::xlim(xlimvec) +
            gtypelabs + ggplot2::theme_classic()

       if (plot) {
            warning("You chose to plot the auc, the figure is not saved.") # nolint
            print(g)
       } else {
            outfile <- paste0("histo_", plottype)
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}
