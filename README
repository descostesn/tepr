# TepR: Transcription elongation profile in R

TepR (Transcription Elongation Profile in R) is an R package designed for analyzing data from nascent RNA sequencing technologies, such as TT-seq, mNET-seq, and PRO-seq.  It calculates the probability distribution of nascent RNA sequencing signal across the gene body or transcription unit of a given gene.  By comparing this profile to a uniform signal, TepR can identify transcription attenuation sites.  Furthermore, it can detect increased or decreased transcription attenuation by comparing profiles across different conditions.  Beyond its rigorous statistical testing and high sensitivity, a key strength of TepR is its ability to resolve the elongation pattern of individual genes, including the precise location of the primary attenuation point, when present.  This capability allows users to visualize and refine genome-wide aggregated analyses, enabling the robust identification of effects specific to gene subsets.  These metrics facilitate comparisons between genes within a condition, across conditions for the same gene, or against a theoretical model of perfect uniform elongation.

## Installation

!!!!!!!!!!!!!!!!!!

## General principles

The preprocessing pipeline consists of the following steps:

1.  Filtering Gencode annotations to extract "transcript" annotations.
2.  Distinguishing between protein-coding (MANE_Select) and long non-coding (lncRNA, Ensembl_canonical) transcripts.
3.  Dividing transcripts into windows of a user-defined size (`windsize`).
4.  Processing bedgraph files to retrieve signal values, excluding blacklisted regions, and retaining scores within high-mappability intervals.
5.  Generating a final annotated table incorporating the scores derived from the preceding steps.

The downstream analysis uses the previously generated final table to:

1. Calculate the average expression levels for transcripts and filter out those below an expression threshold.
2. Evaluate for each transcript the number of sites excluded by the black list and low mappability track.
!!!!!!!!!!!!!


#'   \item countna It takes the result of 'averageandfilterexprs' as input and counts the number of NA values for each transcript based on strand and condition. NA represent missing scores that were filtered out from the black list and mappability track. It returns a data frame where each row corresponds to a transcript, along with its associated gene, strand, and the count of NA values.
#'   \item genesECDF It takes the result of 'averageandfilterexprs' as input and a) filters the main expression table to retain only the expressed transcripts; b) Splits the data by each transcript; c) For each transcript, computes ECDF values for the score columns while respecting the strand orientation ("plus" or "minus"); d) Combines the ECDF results into a final data frame. It returns a list containing two elements: A data frame with ECDF results for each transcript (concatdf) and an integer indicating the number of rows in each transcript table.
#'   \item meandifference It takes the result of 'genesECDF' as input and calculates the mean values, mean Fx (ECDF) and ECDF differences (Fx) for expression data, across different experimental conditions. It returns a data frame that contains, for each condition: mean values for the "value" and "Fx" columns (e.g., \code{mean_value_ctrl}, \code{mean_Fx_ctrl}) and the differences between the \code{Fx} column and coordinate ratios (e.g., \code{diff_Fx_ctrl}).
#'   \item Splitting Split the results of 'meandifference' by transcripts and stores the list into bytranslistmean.
#'   \item allauc It uses the previously computed variable 'bytranslistmean' and it computes the Area Under Curve (AUC) and the differences of AUC between two conditions for a list of transcript data. It returns a data frame containing the AUC and dAUC results for each transcript, along with associated statistical information.
#'   \item kneeid It uses the previously computed variable 'bytranslistmean' and identifies the knee point (i.e., point of maximum change) and the maximum difference in the empirical cumulative distribution function (ECDF) for each transcript, across different experimental conditions. It returns a data frame where each row corresponds to a transcript and contains the coordinates of the knee point and the maximum ECDF difference for each condition.
#'   \item attenuation It uses the results of the previous functions and it computes the attenuation values for each window of each transcript. It returns a data frame containing the computed attenuation values along with associated transcript information.
#'   \item universegroup Using the table produced by 'attenuation', it categorizes genes into a "Universe" and assigns them into groups such as "Attenuated" or "Outgroup" based on transcription data and thresholds. The universe is defined by thresholds for window size, missing data count, mean transcription levels, and p-values. Genes are further classified into groups based on conditions related to AUC and p-value thresholds. A transcript belongs to "Universe" if (window_size > windsizethres & Count_NA < countnathres & meanctrl > meanctrlthres & meanstress > meanstressthres & pvaltheory > pvaltheorythres). A transcript belongs to the groups: - \strong{Attenuated}: if Universe == TRUE & aucstress > aucstressthres & -log10(pvalks) > attenuatedpvalksthres. - \strong{Outgroup}: if Universe == TRUE & pvalks > outgrouppvalksthres & aucctrl > aucctrlthreshigher & aucctrl < aucctrlthreslower.
