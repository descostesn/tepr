.normalizeandsummarize <- function(transvec, dfmeandiff, unigroupdf, daucname, # nolint
    auc_ctrlname, auc_stressname) {

        ## Declaration to tackle CMD check
        transcript <- gene <- strand <- window_size <- coord <- NULL

        ## Selecting full mean and AUC columns
        AUC_allcondi <- unigroupdf %>% dplyr::select(transcript, gene, strand,
            dplyr::contains("Full"), !!sym(daucname), !!sym(auc_ctrlname),
            !!sym(auc_stressname), -contains(c("UP", "DOWN")), window_size)

        ## Selecting expressed transcripts
        res <- dfmeandiff[dfmeandiff$transcript %in% transvec, ]
        ## Join dfmeandiff and AUC_allcondi
        res <- dplyr::left_join(res, AUC_allcondi, by = c("transcript", "gene",
            "strand"))
        ## Selecting columns of interest
        res <- dplyr::select(res, coord, dplyr::contains("mean_value"))
        ## Computing mean for each coordinate
        res <- dplyr::summarise(dplyr::group_by(res, coord), dplyr::across(
            dplyr::contains("mean_value"), \(x) mean(x, na.rm = TRUE)))

        return(res)
}

.checkmetagenes <- function(plottype) {

    if (!isTRUE(all.equal(plottype, "attenuation")) &&
        !isTRUE(all.equal(plottype, "outgroup")) &&
        !isTRUE(all.equal(plottype, "universe")) &&
        !isTRUE(all.equal(plottype, "all")))
        stop("\n[tepr] Error: Invalid plottype.\n",
            "  Use 'attenuation', 'outgroup', 'universe', or 'all'.\n")
}

#' Plot Metagenes for Gene Groups
#'
#' @description
#' This function plots metagene profiles based on transcript data, comparing
#' transcription density across conditions (e.g., control vs. stress). The
#' function allows the user to plot metagenes for different gene groups such as
#' attenuated genes, outgroup genes, the entire universe of genes, or all genes.
#'
#' @usage
#' plotmetagenes(unigroupdf, dfmeandiff, expdf, plottype = "attenuation",
#' daucname = "dAUC_Diff_meanFx_HS_ctrl", auc_ctrlname = "AUC_ctrl",
#' auc_stressname = "AUC_HS", plot = FALSE, formatname = "pdf",
#' outfold = tempdir(), verbose = TRUE)
#'
#' @param unigroupdf A data frame containing gene-level information, including
#'  group classifications and dAUC data for different conditions (see
#'  universegroup).
#' @param dfmeandiff A data frame containing mean transcription values and
#'  coordinates for each transcript (see meandifference).
#' @param expdf A data frame containing experiment data that should have
#'  columns named 'condition', 'replicate', 'strand', and 'path'.
#' @param plottype A string specifying the group of genes to plot. Options are
#'  \code{"attenuation"}, \code{"outgroup"}, \code{"universe"}, or \code{"all"}.
#'  Default is \code{"attenuation"}.
#' @param daucname A string specifying the column name for the delta AUC value
#'  (difference between conditions). Default is
#'  \code{"dAUC_Diff_meanFx_HS_ctrl"}.
#' @param auc_ctrlname A string specifying the column name for the control
#'  condition AUC values. Default is \code{"AUC_ctrl"}.
#' @param auc_stressname A string specifying the column name for the stress
#'  condition AUC values. Default is \code{"AUC_HS"}.
#' @param plot A logical flag indicating whether to display the plot
#'  interactively (\code{TRUE}) or save it to a file (\code{FALSE}).
#'  Default is \code{FALSE}.
#' @param formatname A string specifying the format of the saved plot file.
#'  Default is \code{"pdf"}.
#' @param outfold A string specifying the output folder where the plot will be
#'  saved if \code{plot = FALSE}. Default is the current directory.
#' @param verbose A logical flag indicating whether to display detailed
#'  messages about the function's progress. Default is \code{TRUE}.
#'
#' @return A metagene plot comparing transcription density across conditions
#'  (e.g., control vs. stress) for the selected group of genes. The plot can
#'  either be displayed interactively or saved to a file.
#'
#' @details
#' This function summarizes mean transcription levels across genomic coordinates
#' for different gene groups and plots the transcription density from the
#' transcription start site (TSS) to the transcription termination site (TTS).
#' The function can generate metagene plots for different gene groups such as
#' attenuated, outgroup, or all genes, and compares transcription profiles
#' between conditions (e.g., control vs. stress). The resulting plot helps
#' visualize differences in transcriptional response between groups of genes
#' under different conditions.
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
#'         showtime = FALSE, verbose = FALSE)
#' rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
#' ecdf <- genesECDF(avfilt, verbose = FALSE)
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
#' ## Testing plotmetagenes
#' plotmetagenes(resug, resmeandiff, expdf, plottype = "attenuation", plot = TRUE)
#'
#' @seealso
#' [universegroup], [meandifference]
#'
#' @importFrom ggplot2 ggplot aes geom_line theme_bw ylim labs theme ggsave
#' @importFrom dplyr filter select left_join group_by summarise contains across
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export

plotmetagenes <- function(unigroupdf, dfmeandiff, expdf, plottype = "attenuation",
    daucname = "dAUC_Diff_meanFx_HS_ctrl", auc_ctrlname = "AUC_ctrl",
    auc_stressname = "AUC_HS", plot = FALSE, formatname = "pdf",
    outfold = tempdir(), verbose = TRUE) {

        nbcond <- length(unique(expdf$condition))
        if (!isTRUE(all.equal(nbcond, 2)))
            stop("\n[tepr] Error: Wrong number of conditions.\n",
                "  plotmetagenes requires 2 conditions (found: ", nbcond,
                ").\n")

        .checkmetagenes(plottype)
        colnamevec <- c(daucname, auc_ctrlname, auc_stressname)
        .colnamecheck(colnamevec, unigroupdf)

        if (!file.exists(outfold))
            dir.create(outfold, recursive = TRUE)

        ## Selection of transcripts and define plot title
        if (isTRUE(all.equal(plottype, "attenuation"))) {
            idx <- which(unigroupdf$Group == "Attenuated")
            titleplot <- "Attenuated genes"
        } else if (isTRUE(all.equal(plottype, "outgroup"))) {
            idx <- which(unigroupdf$Group == "Outgroup")
            titleplot <- "Outgroup genes"
        } else if (isTRUE(all.equal(plottype, "universe"))) {
            idx <- which(unigroupdf$Universe)
            titleplot <- "Universe genes"
        } else {
            idx <- seq_len(nrow(unigroupdf))
            titleplot <- "All genes"
        }

        if (isTRUE(all.equal(length(idx), 0)))
            stop("\n[tepr] Error: No transcripts found.\n",
                "  No transcripts match criteria '", plottype, "'.\n")

        transvec <- unigroupdf[idx, "transcript"]
        df <- .normalizeandsummarize(transvec, dfmeandiff, unigroupdf, daucname,
            auc_ctrlname, auc_stressname)
        meanvalctrl <-  colnames(df)[2]
        meanvalstress <- colnames(df)[3]

        ## plotting
        g <-  ggplot2::ggplot() +
            ggplot2::geom_line(data = df, ggplot2::aes(x = .data$coord / 2,
            y = !!sym(meanvalctrl)), color = "#00AFBB", linewidth = 1.5) +
            ggplot2::geom_line(data = df,
                aes(x = .data$coord / 2, y = !!sym(meanvalstress)),
                color = "#FC4E07", linewidth = 1.5) +
            ggplot2::theme_bw() + ggplot2::ylim(0,7) +
            ggplot2::labs(x = "TSS to TTS", title = titleplot,
                subtitle = length(transvec), y = "Transcription density") +
            ggplot2::theme(legend.position = "none", legend.box = "vertical")

        if (plot) {
            warning("[tepr] Warning: Plot displayed only, not saved to file.")
            print(g)
        } else {
            outfile <- paste0("metagene_", plottype)
            if (verbose) message("\t\t Saving plot to ", file.path(outfold,
                paste0(outfile, ".", formatname)))
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                    plot = g, device = formatname, path = outfold)
            }
}
