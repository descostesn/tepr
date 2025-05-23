% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{preprocessing}
\alias{preprocessing}
\title{Preprocess Experimental Data for Genomic Analysis}
\usage{
preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
blacklistpath, genomename = NA, nbcputrans = 1, finaltabpath = tempdir(),
finaltabname = "anno.tsv", tmpfold = file.path(tempdir(), "tmptepr"),
saveobjectpath = tempdir(), savefinaltable = TRUE, reload = FALSE,
showtime = FALSE, showmemory = FALSE, deletetmp = TRUE, chromtab = NA,
forcechrom = FALSE, verbose = TRUE)
}
\arguments{
\item{exptabpath}{Character. Path to the experiment table file.}

\item{gencodepath}{Character. Path to the Gencode annotation file.}

\item{windsize}{Integer. Window size for splitting transcripts.}

\item{maptrackpath}{Character. Path to the mappability track file.}

\item{blacklistpath}{Character. Path to the blacklist file.}

\item{genomename}{Character. Name of the genome assembly (e.g., "hg38").
Default is \code{NA}. If left to NA, chromtab should be provided.}

\item{nbcputrans}{Integer. Number of CPUs to use for transcript processing.
Default is \code{1}.}

\item{finaltabpath}{Character. Path where the final annotated table will be
saved. Default is \code{tempdir()}.}

\item{finaltabname}{Character. Name of the final annotated table file.
Default is \code{anno.tsv}.}

\item{tmpfold}{Character. Path to a temporary folder for intermediate files.
Default is \code{file.path(tempdir(), "tmptepr")}.}

\item{saveobjectpath}{Character. Path to save intermediate objects. Default
is \code{tempdir()}.}

\item{savefinaltable}{Logical. Whether to save the final table to disk.
Default is \code{TRUE}.}

\item{reload}{Logical. Whether to reload intermediate objects if available.
Default is \code{FALSE}.}

\item{showtime}{Logical. Whether to display timing information. Default is
\code{FALSE}.}

\item{showmemory}{Logical. Whether to display memory usage information.
Default is \code{FALSE}.}

\item{deletetmp}{Logical. Whether to delete temporary files after processing.
Default is \code{TRUE}.}

\item{chromtab}{A Seqinfo object retrieved with the rtracklayer method
\code{SeqinfoForUCSCGenome}. If NA, the method is called automatically and
the \code{genomename} should be provided. Default is \code{NA}.}

\item{forcechrom}{Logical indicating if the presence of non-canonical
chromosomes in chromtab (if not NA) should trigger an error. Default is
\code{FALSE}.}

\item{verbose}{Logical. Whether to display detailed progress messages.
Default is \code{TRUE}.}
}
\value{
A data frame representing the final table containing transcript
information and scores on 'windsize' windows for all experiments defined in
the experiment table (exptabpath).
}
\description{
This function orchestrates a pipeline for preprocessing genomic data,
including filtering annotations, splitting transcripts into windows,
retrieving bedgraph values, and generating a final annotated table.
}
\details{
The `preprocessing` function performs several key tasks:
1. Filters Gencode annotations to retrieve "transcript" annotations.
2. Differentiates between protein-coding (MANE_Select) and long non-coding
 (lncRNA, Ensembl_canonical) transcripts.
3. Splits transcripts into windows of size `windsize`.
4. Processes bedgraph files to retrieve values, exclude blacklisted regions,
 and retain high-mappability intervals.
5. Generates a final annotated table with scores derived from the above
steps.

Temporary files created during processing are optionally deleted at the end.
}
\examples{
\donttest{
## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package = "tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
maptrackpath <- system.file("extdata", "k50.umap.chr13.hg38.0.8.bed",
  package = "tepr")
blacklistpath <- system.file("extdata", "hg38-blacklist-chr13.v2.bed",
    package = "tepr")
windsize <- 200
genomename <- "hg38"

## Copying bedgraphs to the current directory
expdfpre <- read.csv(exptabpath)
bgpathvec <- sapply(expdfpre$path, function(x) system.file("extdata", x,
    package = "tepr"))
expdfpre$path <- bgpathvec
write.csv(expdfpre, file = "exptab-preprocessing.csv", row.names = FALSE,
    quote = FALSE)
exptabpath <- "exptab-preprocessing.csv"

## Testing preprocessing
finaltabtest <- preprocessing(exptabpath, gencodepath, windsize, maptrackpath,
    blacklistpath, genomename = genomename, verbose = FALSE)}

}
\seealso{
[retrieveanno], [makewindows], [blacklisthighmap], [createtablescores]
}
