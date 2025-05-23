% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plothistoknee.R
\name{plothistoknee}
\alias{plothistoknee}
\title{Plot Histogram of Distance from TSS to Knee Point}
\usage{
plothistoknee(unigroupdf, plottype = "percent", xlimvec = NA,
binwidthval = NA, kneename = "knee_AUC_HS", plot = FALSE, outfold = tempdir(),
formatname = "pdf", universename = "Universe", groupname = "Group",
verbose = TRUE)
}
\arguments{
\item{unigroupdf}{A data frame containing gene-level statistics, including
knee point data and group classification (see universegroup).}

\item{plottype}{A string specifying the type of distance to plot. Options are
\code{"percent"} for the percentage of the gene or \code{"kb"} for distance
in kilobases. Default is \code{"percent"}.}

\item{xlimvec}{A numeric vector of length 2 specifying the limits of the
x-axis. Default is \code{NA}, which automatically sets the limits based on
\code{plottype}.}

\item{binwidthval}{A numeric value for the width of the bins in the
histogram. Default is \code{NA}, which automatically selects a bin width
based on \code{plottype}.}

\item{kneename}{A string specifying the name of the column in
\code{unigroupdf} that contains the knee point data. Default is
\code{"knee_AUC_HS"}.}

\item{plot}{A logical flag indicating whether to display the plot
interactively (\code{TRUE}) or save it to a file (\code{FALSE}). Default is
\code{FALSE}.}

\item{outfold}{A string specifying the output folder where the plot will be
saved if \code{plot = FALSE}. Default is the current directory.}

\item{formatname}{A string specifying the format of the saved plot file.
Default is \code{"pdf"}.}

\item{universename}{A string specifying the name of the column in
\code{unigroupdf} that defines the universe of genes. Default is
\code{"Universe"}.}

\item{groupname}{A string specifying the name of the column in
\code{unigroupdf} that defines the group classification of genes. Default is
\code{"Group"}.}

\item{verbose}{A logical flag indicating whether to display detailed
messages about the function's progress. Default is \code{TRUE}.}
}
\value{
A histogram showing the distribution of the distance from the TSS to
 the knee point for attenuated genes. The plot can either be displayed
 interactively or saved to a file.
}
\description{
This function generates a histogram showing the distribution of the distance
from the transcription start site (TSS) to the knee point for attenuated
genes. The distance can be plotted either as a percentage of the gene length
or in kilobases (kb).
}
\examples{
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath <- system.file("extdata", "cugusi_6.tsv", package="tepr")
expthres <- 0.1

## Calculating necessary results
expdf <- read.csv(exppath)
transdf <- read.delim(transpath, header = FALSE)
avfilt <- averageandfilterexprs(expdf, transdf, expthres,
        showtime = FALSE, verbose = FALSE)
rescountna <- countna(avfilt, expdf, nbcpu = 1, verbose = FALSE)
ecdf <- genesECDF(avfilt, expdf, verbose = FALSE)
resecdf <- ecdf[[1]]
nbwindows <- ecdf[[2]]
resmeandiff <- meandifference(resecdf, expdf, nbwindows,
    verbose = FALSE)
bytranslistmean <- split(resmeandiff, factor(resmeandiff$transcript))
resknee <- kneeid(bytranslistmean, expdf, verbose = FALSE)
resauc <- allauc(bytranslistmean, expdf, nbwindows, verbose = FALSE)
resatt <- attenuation(resauc, resknee, rescountna, bytranslistmean, expdf,
        resmeandiff, verbose = FALSE)
resug <- universegroup(resatt, expdf, verbose = FALSE)

## Testing plothistoknee
plothistoknee(resug, plot = TRUE)

}
\seealso{
[universegroup]
}
