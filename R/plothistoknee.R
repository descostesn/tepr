#' Plot Histogram of Distance from TSS to Knee Point
#'
#' @description
#' This function generates a histogram showing the distribution of the distance
#' from the transcription start site (TSS) to the knee point for attenuated
#' genes. The distance can be plotted either as a percentage of the gene length
#' or in kilobases (kb).
#'
#' @usage
#' plothistoknee(unigroupdf, plottype = "percent", xlimvec = NA,
#' binwidthval = NA, kneename = "knee_AUC_HS", plot = FALSE, outfold = ".",
#' formatname = "pdf", universename = "Universe", groupname = "Group",
#' verbose = TRUE)
#'
#' @param unigroupdf A data frame containing gene-level statistics, including
#'  knee point data and group classification (see universegroup).
#' @param plottype A string specifying the type of distance to plot. Options are
#'  \code{"percent"} for the percentage of the gene or \code{"kb"} for distance
#'  in kilobases. Default is \code{"percent"}.
#' @param xlimvec A numeric vector of length 2 specifying the limits of the
#'  x-axis. Default is \code{NA}, which automatically sets the limits based on
#'  \code{plottype}.
#' @param binwidthval A numeric value for the width of the bins in the
#'  histogram. Default is \code{NA}, which automatically selects a bin width
#'  based on \code{plottype}.
#' @param kneename A string specifying the name of the column in
#'  \code{unigroupdf} that contains the knee point data. Default is
#'  \code{"knee_AUC_HS"}.
#' @param plot A logical flag indicating whether to display the plot
#'  interactively (\code{TRUE}) or save it to a file (\code{FALSE}). Default is
#'  \code{FALSE}.
#' @param outfold A string specifying the output folder where the plot will be
#'  saved if \code{plot = FALSE}. Default is the current directory.
#' @param formatname A string specifying the format of the saved plot file.
#'  Default is \code{"pdf"}.
#' @param universename A string specifying the name of the column in
#'  \code{unigroupdf} that defines the universe of genes. Default is
#'  \code{"Universe"}.
#' @param groupname A string specifying the name of the column in
#'  \code{unigroupdf} that defines the group classification of genes. Default is
#'  \code{"Group"}.
#' @param verbose A logical flag indicating whether to display detailed
#'  messages about the function's progress. Default is \code{TRUE}.
#'
#' @return A histogram showing the distribution of the distance from the TSS to
#'  the knee point for attenuated genes. The plot can either be displayed
#'  interactively or saved to a file.
#'
#' @examples
#' # Assuming `unigroupdf` contains the necessary data:
#' # plothistoknee(unigroupdf, plottype = "kb", xlimvec = c(0, 300),
#' # binwidthval = 10)
#'
#' @seealso
#' [universegroup]
#'
#' @importFrom ggplot2 ggplot aes geom_histogram xlab xlim labs theme_classic ggsave
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export

plothistoknee <- function(unigroupdf, plottype = "percent", xlimvec = NA, # nolint
    binwidthval = NA, kneename = "knee_AUC_HS", plot = FALSE, outfold = ".",
    formatname = "pdf", universename = "Universe", groupname = "Group",
    verbose = TRUE) {

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
            gtypeaes <- ggplot2::aes(x = (!!sym(kneename) *
                .data$window_size) / 1000)
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
            outfile <- paste0("histo_", plottype, ".", formatname)
            if (verbose) message("Saving plot to ", file.path(outfold, outfile))
            ggplot2::ggsave(filename = outfile, plot = g, device = formatname,
                path = outfold)
        }
}
