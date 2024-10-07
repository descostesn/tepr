.normalizeandsummarize <- function(transvec, dfmeandiff, unigroupdf, daucname, # nolint
    auc_ctrlname, auc_stressname) {

    ## Selecting full mean and AUC columns
    AUC_allcondi <- unigroupdf %>% dplyr::select(transcript, gene, strand,
        dplyr::contains("Full"), !!sym(daucname), !!sym(auc_ctrlname),
        !!sym(auc_stressname), -contains(c("UP", "DOWN")), window_size)

    ## Selecting coord and mean values
    result <- dfmeandiff %>%
        dplyr::filter(transcript %in% transvec) %>% #nolint
        dplyr::left_join(., AUC_allcondi,
            by = c("transcript", "gene")) %>%
        dplyr::select(transcript, gene, coord, dplyr::contains("mean_value"),
        -dplyr::contains("Full"))  %>% dplyr::group_by(coord) %>%
        dplyr::summarise(dplyr::across(dplyr::contains("mean_value"),
        ~ mean(., na.rm = TRUE)))

    return(result)
}

.checkmetagenes <- function(plottype) {

    if (!isTRUE(all.equal(plottype, "attenuation")) &&
        !isTRUE(all.equal(plottype, "outgroup")) &&
        !isTRUE(all.equal(plottype, "universe")) &&
        !isTRUE(all.equal(plottype, "all")))
        stop("plot type should be one of: attenuation, outgroup, universe, all")
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
#' plotmetagenes(unigroupdf, dfmeandiff, plottype = "attenuation",
#' daucname = "dAUC_Diff_meanFx_HS_ctrl", auc_ctrlname = "AUC_ctrl",
#' auc_stressname = "AUC_HS", plot = FALSE, formatname = "pdf", outfold = ".",
#' verbose = TRUE)
#'
#' @param unigroupdf A data frame containing gene-level information, including
#'  group classifications and dAUC data for different conditions (see
#'  universegroup).
#' @param dfmeandiff A data frame containing mean transcription values and
#'  coordinates for each transcript (see meandifference).
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
#' # Assuming `unigroupdf` and `dfmeandiff` contain the necessary data:
#' # plotmetagenes(unigroupdf, dfmeandiff, plottype = "universe", plot = TRUE)
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

plotmetagenes <- function(unigroupdf, dfmeandiff, plottype = "attenuation",
    daucname = "dAUC_Diff_meanFx_HS_ctrl", auc_ctrlname = "AUC_ctrl",
    auc_stressname = "AUC_HS", plot = FALSE, formatname = "pdf",
    outfold = ".", verbose = TRUE) {

    .checkmetagenes(plottype)
    colnamevec <- c(daucname, auc_ctrlname, auc_stressname)
    .colnamecheck(colnamevec, unigroupdf)

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
        stop("No transcripts were found for the criteria ", plottype)

     transvec <- unigroupdf[idx, "transcript"]
     df <- .normalizeandsummarize(transvec, dfmeandiff, unigroupdf, daucname,
        auc_ctrlname, auc_stressname)
    meanvalctrl <-  colnames(df)[2]
    meanvalstress <- colnames(df)[3]

    ## plotting
    g <-  ggplot2::ggplot() +
        ggplot2::geom_line(data = df, ggplot2::aes(x = .data$coord / 2,
        y = !!sym(meanvalctrl)), color = "#00AFBB", size = 1.5) +
        ggplot2::geom_line(data = df,
            aes(x = .data$coord / 2, y = !!sym(meanvalstress)),
            color = "#FC4E07", size = 1.5) +
        ggplot2::theme_bw() + ggplot2::ylim(0,7) +
        ggplot2::labs(x = "TSS to TTS", title = titleplot,
            subtitle = length(transvec), y = "Transcription density") +
        ggplot2::theme(legend.position = "none", legend.box = "vertical")

    if (plot) {
        warning("You chose to plot the auc, the figure is not saved.") # nolint
        print(g)
    } else {
        outfile <- paste0("metagene_", plottype)
        if (verbose) message("\t\t Saving plot to ", file.path(outfold,
            paste0(outfile, ".", formatname)))
        ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}
