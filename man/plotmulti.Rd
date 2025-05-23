% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotmulti.R
\name{plotmulti}
\alias{plotmulti}
\title{Generate all tepr plots for all experiment comparisons}
\usage{
plotmulti(resteprmulti, expdf, ecdfgenevec, outfold = tempdir(), digits = 2,
middlewind = 100, pval = 0.01, colvec = c("#90AFBB", "#10AFBB",
"#FF9A04", "#FC4E07"), genaucvec = NA, aucaxisminx = -10,
aucaxismaxx = 100, aucaxisminy = -10, aucaxismaxy = 100, aucmaintitle = "",
aucsubtitle = "", auclegendpos = "bottom", formatname = "pdf",
uniname = "Universe", groupname = "Group", histkneexlim = NA,
binwidthvalhistknee = NA, verbose = TRUE)
}
\arguments{
\item{resteprmulti}{Result returned by the function \code{teprmulti}.}

\item{expdf}{A data frame containing experiment data that should have
columns named 'condition', 'replicate', 'strand', and 'path'.}

\item{ecdfgenevec}{A vector specifying the names of the genes of interest to
plot the ecdf of.}

\item{outfold}{Path to the output folder where the plots will be written.
Subfolders with the names of the comparisons are automatically created.
Default is \code{tempdir()}.}

\item{digits}{For the ecdf plot, the number of decimal places to round the
AUC and KS values. Default is \code{2}.}

\item{middlewind}{For the ecdf plot, the index of the middle window
representing the region centered around the TSS. Default is \code{100}.}

\item{pval}{For the ecdf plot, a numeric value for the p-value threshold to
determine the significance of the KS test. Default is \code{0.01}.}

\item{colvec}{For the ecdf plot, a vector of 4 colors used to distinguish
the different conditions. Default is \code{c("#90AFBB", "#10AFBB",
"#FF9A04", "#FC4E07")}.}

\item{genaucvec}{For the auc plot, vector of gene names to highlight,
Used for the plot of type "pval". Default is \code{NA}. If left to NA, the
plot taking into account the p-values is not generated.}

\item{aucaxisminx}{For the auc plot, minimum value for the x-axis. Default
is \code{-10}.}

\item{aucaxismaxx}{For the auc plot, maximum value for the x-axis. Default
is \code{100}.}

\item{aucaxisminy}{For the auc plot, minimum value for the y-axis. Default
is \code{-10}.}

\item{aucaxismaxy}{For the auc plot, maximum value for the y-axis. Default
is \code{100}.}

\item{aucmaintitle}{For the auc plot, main title of the plot. Default is an
empty string.}

\item{aucsubtitle}{For the auc plot, subtitle of the plot. Default is an
empty string.}

\item{auclegendpos}{For the auc plot, position of the legend. Default is
\code{"bottom"}.}

\item{formatname}{Format of the saved plot (e.g., "pdf", "png"). Default is
\code{"pdf"}.}

\item{uniname}{Column name in the second element of \code{resteprmulti}
representing the universe selection. Default is \code{"Universe"}.}

\item{groupname}{Column name in the second element of \code{resteprmulti}
 representing the type of group a transcript belong to. Default is
\code{"Group"}.}

\item{histkneexlim}{For the plot histoknee, a numeric vector of length 2
 specifying the limits of the x-axis. Default is \code{NA}, which
automatically sets the limits based on \code{plottype}.}

\item{binwidthvalhistknee}{For the plot histoknee, a numeric value for the
width of the bins in the histogram. Default is \code{NA}, which
automatically selects a bin width based on \code{plottype}.}

\item{verbose}{A logical flag indicating whether to display detailed
messages about the function's progress. Default is \code{TRUE}.}
}
\value{
Nothing is returned. Figures are written to outfold in the subfolder
 of the corresponding comparison.
}
\description{
This function generates for all experiment comparisons contained in the
object \code{resteprmulti} all plots of tepr: ECDF, auc, metagene, and
histtoknee.
}
\details{
The function goes through each element of resteprmulti which corresponds to
a comparison of two conditions. For each element it calls the following
functions:
#' \itemize{
  \item \code{"plotecdf"}: The function generates a figure for each gene
 given in ecdfgenevec.
  \item \code{"plotauc"}: Generates figures by groups and pval. The lattest
 figure is not generated if genaucvec = NA.
  \item \code{"plotmetagenes"}: Generates the figures by \code{attenuation},
 \code{outgroup}, \code{universe}, and \code{all}.
  \item \code{"plothistoknee"}: Generate the figures by \code{percent} and
 \code{kb}.
}
}
\examples{
# Assuming resteprmulti is the object returned by the function teprmulti
# and expdf contains the necessary data:
\dontrun{
plotmulti(resteprmulti, expdf, ecdfgenevec = c("EGFR", "DAP", "FLI1"))}

}
\seealso{
[teprmulti], [plotecdf], [plotauc], [plotmetagenes], [plothistoknee]
}
