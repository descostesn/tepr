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
    res <- list(subtext, vlinedf, kneeval)
    return(res)
}

.valcolbuild <- function(condvec, repvec) {
    return(as.vector(sapply(condvec,
        function(cond) {
            sapply(repvec, function(rep) {
            paste("value", cond, paste0("rep", rep), "score", sep = "_") })})))
}

.callggplotecdf <- function(dflongecdf, colvec, windsizefact, vlinedf, subtext,
    outfold, genename, kneeval, plot, formatname, verbose) {

    colvec <- as.vector(factor(dflongecdf$conditions, labels = colvec))
    ylimval <- 2 * max(dflongecdf$value)
    linexvals <- dflongecdf$coord * windsizefact
    lineyvals <- dflongecdf$coord / max(dflongecdf$coord)
    areayvals <- dflongecdf$value / ylimval

    g <- ggplot2::ggplot(dflongecdf, ggplot2::aes(x = .data$coord, y = .data$Fx,
        color = .data$conditions)) + # nolint
        ggplot2::geom_line(aes(x = linexvals, y = lineyvals),
            linetype = "dashed", color = "red")

    g1 <- g + ggplot2::geom_area(ggplot2::aes(x = linexvals, y = areayvals,
        fill = .data$conditions), # nolint
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
        if (verbose) message("\t\t Saving figure to ", file.path(outfold,
            paste0(genename, ".", formatname)))
        ggplot2::ggsave(filename = paste0(genename, ".", formatname),
            plot = g2, device = formatname, path = outfold)
    }
}

.windowsizefactor <- function(df, middlewind) {

    idxmiddle <- which(df$coord == middlewind)
    if (isTRUE(all.equal(length(idxmiddle), 0)))
        stop("\n\t The window ", middlewind, " was not found. Did you define a",
        " number of windows equal to ", middlewind*2, "? If not, adjust the ",
        "parameter 'middlewind' to the half of the number of windows.\n")
    dfmid <- df[idxmiddle, ]
    windsizefact <- (dfmid$coor2 - dfmid$coor1) / 1000
    return(windsizefact)
}

#' Plot Empirical Cumulative Distribution Function (ECDF)
#'
#' @description
#' This function generates an ECDF plot to analyze transcription density
#' relative to the distance from the transcription start site (TSS) across
#' different conditions. The plot displays AUC values, Kolmogorov-Smirnov (KS)
#' statistics, and knee points, with options to display or save the plot.
#'
#' @usage
#' plotecdf(dfmeandiff, unigroupdf, expdf, genename,
#'    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
#'    outfold = tempdir(), digits = 2, middlewind = 100, pval = 0.01,
#'    plot = FALSE, formatname = "pdf", verbose = TRUE)
#'
#'
#' @param dfmeandiff A data frame containing the mean differences of
#'  transcription levels and cumulative distribution values (Fx) for different
#'  windows around the TSS (see meandifference).
#' @param unigroupdf A data frame containing gene-specific statistics, including
#'  their belonging to Universe or Group (see universegroup).
#' @param expdf A data frame containing experiment data that should have
#'              columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param genename A string specifying the name of the gene of interest to plot.
#' @param colvec A vector of colors used to distinguish different conditions in
#'  the plot. Default is \code{c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")}.
#' @param outfold A string specifying the output folder where the plot will be
#'  saved if \code{plot = FALSE}. Default is \code{tempdir()}.
#' @param digits The number of decimal places to round the AUC and KS values.
#'  Default is \code{2}.
#' @param middlewind The index of the middle window representing the region
#'  centered around the TSS. Default is \code{100}.
#' @param pval A numeric value for the p-value threshold to determine the
#'  significance of the KS test. Default is \code{0.01}.
#' @param plot A logical flag indicating whether to display the plot
#'  interactively (\code{TRUE}) or save it to a file (\code{FALSE}). Default is
#'  \code{FALSE}.
#' @param formatname String of the format of the saved plot. Possible values are
#'  "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", and "svg".
#'  Default is \code{"pdf"}.
#' @param verbose A logical flag indicating whether to display detailed
#'  messages about the function's progress. Default is \code{TRUE}.
#'
#' @return An ECDF plot showing the transcription density across windows around
#'  the TSS, with highlights for significant KS test results and knee points.
#' The plot can either be displayed or saved as a file.
#'
#' @details
#' The function processes data related to transcription levels and cumulative
#' transcription density for a given gene across multiple experimental
#' conditions. The ECDF plot is constructed with optional annotation of key
#' statistics such as AUC values and significant KS test results. Knee points,
#' representing significant changes in transcription density, are also displayed
#' if the KS test passes the specified p-value threshold.
#'
#' Colvec: The number of colors should be equal to the number of rows of expdf
#' divided by two (a forward and reverse files are provided for each
#' experiment).
#'
#' @examples
#' exppath <-  system.file("extdata", "exptab.csv", package="tepr")
#' transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
#' expthres <- 0.1
#'
#' ## Calculating necessary results
#' expdf <- read.csv(exppath)
#' transdf <- read.delim(transpath, header = FALSE)
#' avfilt <- averageandfilterexprs(expdf, transdf, expthres,
#'        showtime = FALSE, verbose = FALSE)
#' rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' resmeandiff <- meandifference(resecdf, expdf, nbwindows,
#'     verbose = FALSE)
#' bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
#' resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
#' resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
#' resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
#'         resmeandiff, verbose = FALSE)
#' resug <- universegroup(resatt, expdf, verbose = FALSE)
#'
#' ## Testing plotecdf
#' colvec <- c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
#' plotecdf(resmeandiff, resug, expdf, "EGFR", colvec, plot = TRUE, verbose = FALSE)
#'
#' @seealso
#' [meandifference], [universegroup]
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_area geom_vline scale_fill_manual scale_color_manual labs theme_classic sec_axis ggsave
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

plotecdf <- function(dfmeandiff, unigroupdf, expdf, genename,  # nolint
    colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07"),
    outfold = tempdir(), digits = 2, middlewind = 100, pval = 0.01, plot = FALSE,
    formatname = "pdf", verbose = TRUE) {

        nbrep <- length(expdf$replicate) / 2
        if (!isTRUE(all.equal(length(colvec), nbrep)))
            stop("\n\t The vector of colours colvec should have ", nbrep,
                " values.\n")

        if (verbose) message("\n Plotting ecdf for gene ", genename)

        if (!file.exists(outfold))
            dir.create(outfold, recursive = TRUE)

        ## Retrieving rows concerning the gene of interest
        if (verbose) message("\t Retrieving rows concerning the gene of ",
            "interest")
        idxgene <- which(dfmeandiff$gene == genename)
        if (isTRUE(all.equal(length(idxgene), 0)))
            stop("\n\t The gene ", genename, " was not found.\n")
        df <- dfmeandiff[idxgene, ]
        idxinfo <- which(unigroupdf$gene == genename)
        if (isTRUE(all.equal(length(idxinfo), 0)))
            stop("\n\t The gene ", genename, " was not found in unigroupdf.\n")
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
            dplyr::mutate(conditions = gsub("Fx_|_score", "",
                .data$conditions))
        dflongval <- df %>% tidyr::pivot_longer(
            cols = tidyselect::all_of(valcolvec), names_to = "conditions",
            values_to = "value") %>%
            dplyr::mutate(conditions = gsub("value_|_score", "",
                .data$conditions))
        ## merging
        commoncols <- intersect(names(dflongfx), names(dflongval))
        dflongecdf <- merge(dflongfx, dflongval, by = commoncols)

        ## Plotting
        if (verbose && !plot) message("\t Generating ecdf plot to ", outfold)
        .callggplotecdf(dflongecdf, colvec, windsizefact, vlinedf, subtext,
            outfold, genename, kneeval, plot, formatname, verbose)
}
