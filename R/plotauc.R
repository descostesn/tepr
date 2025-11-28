.callggplotauc <- function(df, aesvar, geompointinfo, geompointinfo2,
    geompointinfo3, plottype, axismin_x, axismax_x, axismin_y, axismax_y,
    labelx, labely, maintitle, subtitle, legendpos, plot, outfile, formatname,
    outfold, genevec, verbose) {

        ## Declaration to tackle CMD check
        gene <- NULL

        ## Structure of the basic scatterplot
        g <- ggplot2::ggplot(df, aesvar) + geompointinfo + geompointinfo2

        if (isTRUE(all.equal(plottype, "pval")) && !is.na(genevec[1])) {

            ## Adding highlight of the genes
            g <- g + ggrepel::geom_label_repel(data = subset(df,
                gene %in% genevec), aes(label = .data$gene), # nolint
                box.padding = 0.55, point.padding = 0,
                segment.color = "black", max.overlaps = 50, color = "red") +
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
            warning("[tepr] Warning: Plot displayed only, not saved to file.")
            print(g)
        } else {
            if (verbose) message("\t\t Saving plot to ", file.path(outfold,
                paste0(outfile, ".", formatname)))
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}

.checkplotaucparams <- function(plottype, auc_ctrlname, auc_stressname,
    pvalkstestcolname, genevec, tab, expdf) {

        nbcond <- length(unique(expdf$condition))
        if (!isTRUE(all.equal(nbcond, 2)))
            stop("\n[tepr] Error: Wrong number of conditions.\n",
                "  plotauc requires exactly 2 conditions (found: ", nbcond,
                ").\n")

        if (!isTRUE(all.equal(plottype, "pval")) &&
            !isTRUE(all.equal(plottype, "groups")))
                stop("\n[tepr] Error: Invalid plottype.\n",
                    "  Use 'pval' or 'groups'.\n")

        colnamevec <- c(auc_ctrlname, auc_stressname, pvalkstestcolname)
        .colnamecheck(colnamevec, tab)

        if (isTRUE(all.equal(plottype, "groups")) && !is.na(genevec[1]))
            stop("\n[tepr] Error: Unnecessary parameter.\n",
                "  'genevec' is not used for plottype='groups'.\n")
}

#' Plot AUC Comparison Between Conditions
#'
#' @description
#' This function generates scatterplots comparing the area under the curve (AUC)
#' for control and stress conditions, with an option to highlight specific genes
#' or groups. The plot can be saved as a file or displayed interactively.
#'
#' @usage
#' plotauc(tab, expdf, genevec = NA, auc_ctrlname = "AUC_ctrl",
#' auc_stressname = "AUC_HS",
#' pvalkstestcolname = "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
#' labelx = "AUC in Control", labely = "AUC in Stress", axismin_x = -10,
#' axismax_x = 100, axismin_y = -10, axismax_y = 100, maintitle = "",
#' subtitle = "", legendpos = "bottom", formatname = "pdf", outfold = tempdir(),
#' outfile = "AUCcompare_pval", plottype = "pval", plot = FALSE,
#' universename = "Universe", groupname = "Group", verbose = TRUE)
#'
#' @param tab A data frame containing the AUC values for control and stress
#'  conditions, and other columns required for plotting (e.g., p-values or
#'  group memberships, see allauc).
#' @param expdf A data frame containing experiment data that should have
#'  columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param genevec A vector of gene names to highlight on the plot, applicable
#'  when \code{plottype} is set to "pval". Default is \code{NA}.
#' @param auc_ctrlname The column name in \code{tab} for the AUC under control
#'  conditions. Default is \code{"AUC_ctrl"}.
#' @param auc_stressname The column name in \code{tab} for the AUC under stress
#'  conditions. Default is \code{"AUC_HS"}.
#' @param pvalkstestcolname The column name in \code{tab} for the adjusted FDR
#'  p-values from the KS test. Default is
#'  \code{"adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"}.
#' @param labelx Label for the x-axis. Default is \code{"AUC in Control"}.
#' @param labely Label for the y-axis. Default is \code{"AUC in Stress"}.
#' @param axismin_x Minimum value for the x-axis. Default is \code{-10}.
#' @param axismax_x Maximum value for the x-axis. Default is \code{100}.
#' @param axismin_y Minimum value for the y-axis. Default is \code{-10}.
#' @param axismax_y Maximum value for the y-axis. Default is \code{100}.
#' @param maintitle Main title of the plot. Default is an empty string.
#' @param subtitle Subtitle of the plot. Default is an empty string.
#' @param legendpos Position of the legend. Default is \code{"bottom"}.
#' @param formatname Format of the saved plot (e.g., "pdf", "png"). Default is
#'  \code{"pdf"}.
#' @param outfold Output folder where the plot will be saved. Default is
#'  \code{tempdir()}.
#' @param outfile Name of the output file. Default is
#'  \code{"AUCcompare_pval"}.
#' @param plottype Type of plot to generate. Can be \code{"pval"} for p-value
#'  based plots or \code{"groups"} for group-based plots. Default is
#'  \code{"pval"}.
#' @param plot A logical flag indicating whether to display the plot
#'  interactively (\code{TRUE}) or save it to a file (\code{FALSE}). Default is
#'  \code{FALSE}.
#' @param universename Column name in \code{tab} representing the universe
#'  group in group-based plots. Default is \code{"Universe"}.
#' @param groupname Column name in \code{tab} representing specific groups in
#'  group-based plots. Default is \code{"Group"}.
#' @param verbose A logical flag indicating whether to display detailed
#'  messages about the function's progress. Default is \code{TRUE}.
#'
#' @return A plot comparing AUC values between control and stress conditions,
#'  either displayed or saved to a file.
#'
#' @details
#' The function supports two plot types:
#' \itemize{
#'   \item \code{"pval"}: The plot highlights genes based on adjusted FDR
#'  p-values and can highlight specific genes provided in \code{genevec}.
#'   \item \code{"groups"}: The plot highlights predefined groups, such as
#'  "Attenuated" and "Outgroup", within the data.
#' }
#'
#' If \code{plot = TRUE}, the plot is displayed interactively. If
#'  \code{plot = FALSE}, the plot is saved to a file in the specified format and
#'  output folder.
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
#' ecdf <- genesECDF(avfilt, verbose = FALSE)
#' resecdf <- ecdf[[1]]
#' nbwindows <- ecdf[[2]]
#' resmeandiff <- meandifference(resecdf, expdf, nbwindows,
#'    verbose = FALSE)
#' bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
#' resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
#' resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
#' resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
#'        resmeandiff, verbose = FALSE)
#' resug <- universegroup(resatt, expdf, verbose = FALSE)
#'
#' ## Testing plotauc
#' plotauc(resug, expdf, plottype = "groups", plot = TRUE)
#'
#' @seealso
#' [allauc]
#'
#' @importFrom dplyr arrange filter
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d labs coord_fixed theme_classic theme xlim ylim ggsave
#' @importFrom ggrepel geom_label_repel
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom MASS kde2d
#'
#' @export

plotauc <- function(tab, expdf, genevec = NA, # nolint
    auc_ctrlname = "AUC_ctrl", auc_stressname = "AUC_HS",
    pvalkstestcolname = "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in Stress", axismin_x = -10,
    axismax_x = 100, axismin_y = -10, axismax_y = 100, maintitle = "",
    subtitle = "", legendpos = "bottom", formatname = "pdf", outfold = tempdir(),
    outfile = "AUCcompare_pval", plottype = "pval", plot = FALSE,
    universename = "Universe", groupname = "Group", verbose = TRUE) {

        .checkplotaucparams(plottype, auc_ctrlname, auc_stressname,
            pvalkstestcolname, genevec, tab, expdf)

        if (!file.exists(outfold))
            dir.create(outfold, recursive = TRUE)

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
            outfile, formatname, outfold, genevec, verbose)
}
