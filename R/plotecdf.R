.subtextvline <- function(condvec, geneinfo, digits, pval) {

    reslist <- lapply(condvec, function(cond, geneinfo, digits) {

        ksname <- paste0("adjFDR_p_AUC_", cond)
        ksval <- round(geneinfo[, ksname], digits)
        aucval <- round(geneinfo[, paste0("AUC_", cond)], digits)

        if (ksval < pval) {
            kneeaucval <- geneinfo[, paste0("knee_AUC_", cond)]
            kneeval <- round(kneeaucval * geneinfo$window_size / 1000,
                digits)
            vlinecond <- data.frame(condition = cond, kneeval = kneeval)
        } else {
            kneeval <- "NA"
            vlinecond <- data.frame(condition = character(),
                kneeval = numeric())
        }
        condtext <- paste0(cond, ": AUC = ", aucval, ", KS = ", ksval,
            ", Knee (kb) = ", kneeval)
        return(list(condtext, vlinecond, kneeval))
    }, geneinfo, digits)
    subtext <- paste(sapply(reslist, "[", 1), collapse = "\n")
    vlinedf <- do.call("rbind", sapply(reslist, "[", 2))
    kneeval <- sapply(reslist, "[", 3)
    return(list(subtext, vlinedf, kneeval))
}

.valcolbuild <- function(condvec, repvec) {
    return(as.vector(sapply(condvec,
        function(cond) {
            sapply(repvec, function(rep) {
            paste("value", cond, paste0("rep", rep), "score", sep = "_") })})))
}

.callggplotecdf <- function(dflongecdf, colvec, windsizefact, vlinedf, subtext,
    outfold, genename, kneeval, plot) {

    colvec <- as.vector(factor(dflongecdf$conditions, labels = colvec))
    ylimval <- 2 * max(dflongecdf$value)
    linexvals <- dflongecdf$coord * windsizefact
    lineyvals <- dflongecdf$coord / max(dflongecdf$coord)
    areayvals <- dflongecdf$value / ylimval

    g <- ggplot2::ggplot(dflongecdf, ggplot2::aes(x = coord, y = Fx, # nolint
        color = conditions)) + # nolint
        ggplot2::geom_line(aes(x = linexvals, y = lineyvals),
            linetype = "dashed", color = "red")

    g1 <- g + ggplot2::geom_area(ggplot2::aes(x = linexvals, y = areayvals,
        fill = conditions), # nolint
        alpha = 0.1, linewidth = 0.2, position = 'identity') + # nolint
        ggplot2::scale_fill_manual(values = colvec) +
        ggplot2::geom_line(linewidth = 1, aes(x = linexvals)) +
        ggplot2::scale_color_manual(values = colvec)

    g2 <- g1 + ggplot2::scale_y_continuous(sec.axis =
        ggplot2::sec_axis(~ . * ylimval, name = "Transcription level")) +
        ggplot2::labs(x = "Distance from TSS (kb)",
             y = "Cumulative transcription density", title = genename,
             subtitle = subtext) + ggplot2::theme_classic()

    if (!isTRUE(all.equal(nrow(vlinedf), 0)))
        g2 <- g2 + ggplot2::geom_vline(data = vlinedf,
            ggplot2::aes(xintercept = kneeval),
            linetype = "dashed", color = "darkgrey")

    if (plot) {
        warning("You chose to plot the ecdf, the figure is not saved.")
        print(g2)
    } else {
        ggplot2::ggsave(filename = paste0(genename, ".pdf"),
            plot = g2, device = "pdf", path = outfold)
    }
}

.windowsizefactor <- function(df, middlewind) {

    idxmiddle <- which(df$coord == middlewind)
    if (isTRUE(all.equal(length(idxmiddle), 0)))
        stop("The window ", middlewind, " was not found. Did you define a ",
        "number of windows equal to ", middlewind*2, "? If not, adjust the ",
        "parameter 'middlewind' to the half of the number of windows.")
    dfmid <- df[idxmiddle, ]
    windsizefact <- (dfmid$coor2 - dfmid$coor1) / 1000
    return(windsizefact)
}

plotecdf <- function(dfmeandiff, unigroupdf, expdf, genename, colvec, outfold, # nolint
    digits = 2, middlewind = 100, pval = 0.01, plot = FALSE, verbose = TRUE) {

    ## Retrieving rows concerning the gene of interest
    if (verbose) message("\t Retrieving rows concerning the gene of interest")
    idxgene <- which(dfmeandiff$gene == genename)
    if (isTRUE(all.equal(length(idxgene), 0)))
        stop("The gene ", genename, " was not found")
    df <- dfmeandiff[idxgene, ]
    idxinfo <- which(unigroupdf$gene == genename)
    if (isTRUE(all.equal(length(idxinfo), 0)))
        stop("The gene ", genename, " was not found in unigroupdf")
    geneinfo <- unigroupdf[idxinfo, ]

    if (verbose) message("\t Gathering statistics about each condition")
    ## Computing the window size factor
    windsizefact <- .windowsizefactor(df, middlewind)

    ## Retrieving auc, ks, and knee
    condvec <- unique(expdf$condition)
    restmp <- .subtextvline(condvec, geneinfo, digits, pval)
    subtext <- restmp[[1]]
    vlinedf <- restmp[[2]]
    kneeval <- restmp[[3]]

    ## Building data.frame for plotting with fx and value
    if (verbose) message("\t Building df for plotting with fx and value")
    repvec <- unique(expdf$replicate)
    colnamedfvec <- colnames(df)
    fxcolvec <- colnamedfvec[grep("^Fx_", colnamedfvec)]
    valcolvec <- .valcolbuild(condvec, repvec)
    ## Apply pivot
    dflongfx <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(fxcolvec), names_to = "conditions",
        values_to = "Fx") %>%
        dplyr::mutate(conditions = gsub("Fx_|_score", "", conditions)) # nolint
    dflongval <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(valcolvec), names_to = "conditions",
        values_to = "value") %>%
        dplyr::mutate(conditions = gsub("value_|_score", "", conditions)) # nolint
    ## merging
    commoncols <- intersect(names(dflongfx), names(dflongval))
    dflongecdf <- merge(dflongfx, dflongval, by = commoncols)

    ## Plotting
    if (verbose && !plot) message("\t Generating ecdf plot to ", outfold)
    .callggplotecdf(dflongecdf, colvec, windsizefact, vlinedf, subtext, outfold,
        genename, kneeval, plot)
}
