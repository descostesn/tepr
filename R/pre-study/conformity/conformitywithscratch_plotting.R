library("ggplot2")
library("tidyr")
library("tidyselect")

## /g/romebioinfo/tmp/comparewithscratch-plotting


##################
# PARAMETERS
##################

unigroupdfpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_unigroupdf.rds" # nolint
dfmeandiffpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_dfmeandiffvic.rds" # nolint
expdfpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
colvec <- c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
outfold <- "/g/romebioinfo/tmp/comparewithscratch-plotting"


##################
#FUNCTIONS
##################

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

        colnamecheckvec <- c(auc_ctrlname, auc_stressname, pvalkstestcolname)
        invisible(sapply(colnamecheckvec, function(currentcol, tab) {
            idx <- grep(currentcol, colnames(tab))
            if (isTRUE(all.equal(length(idx), 0)))
                stop("The column ", currentcol, " does not exist in the ",
                    "provided table.")
        }, tab))

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



##################
# MAIN
##################

## Reading objects and files
dfmeandiff <- readRDS(dfmeandiffpath)
unigroupdf <- readRDS(unigroupdfpath)
expdf <- read.csv(expdfpath, header = TRUE)


####
#### plotecdf
####

plotecdf(dfmeandiff, unigroupdf, expdf, "EGFR", colvec, outfold, plot = TRUE)
plotecdf(dfmeandiff, unigroupdf, expdf, "MARCHF6", colvec, outfold, plot = TRUE)


####
#### plotauc
####

genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")
plotauc(unigroupdf, genevec, plot = TRUE)
plotauc(unigroupdf, legendpos = "none", subtitle = "Genes selected for Unibind",
    maintitle = "AUC Control vs HS", plot = TRUE, plottype = "groups")
