####################################
# This script goes through documentation/explore.R and homogenizes it with
# preprocessing.R
#
# Descostes - R-4.4.1 - July 2024
####################################

library("tidyr")
library("dplyr")
library("tidyselect")

##################
# PARAMETERS
##################

alldfpath <- "/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv"
expthres <- 0.1
nbcpu <- 6


##################
# MAIN
##################

## Reading alldf and info tab
alldf <- readRDS(alldfpath)
expdf <- read.csv(exptabpath, header = TRUE)

## Filtering out non expressed transcripts:
## 1) for each column, calculate the average expression per transcript (over each frame) # nolint
## 2) For each column, remove a line if it contains only values < expthres separating strands # nolint


!!!!!!!!!!!!!!!!!!!!!
scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction)
idxscores <- sapply(scorecolvec, grep, colnames(alldf))

dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
    dplyr::summarize(gene = gene[1], strand = strand[1],
        dplyr::across(
            tidyselect::all_of(scorecolvec),
            ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint

dfstrandlist <- mapply(function(strandname, directname, dfbytrans, expthres) {
    if ((isTRUE(all.equal(strandname, "+")) && isTRUE(all.equal(directname, "rev"))) || 
        (isTRUE(all.equal(strandname, "-")) && isTRUE(all.equal(directname, "fwd"))))
        stop("Strand and direction do not match, contact the developper")

    dfstrand <- dfbytranscript %>%
    dplyr::filter(strand == strandname) %>%
    select(gene, transcript, strand, tidyselect::contains(directname))  %>%
    dplyr::filter(dplyr::if_all(tidyselect::all_of(tidyselect::contains("mean")), ~ !is.na(.))) %>%
    dplyr::filter(dplyr::if_all(tidyselect::all_of(tidyselect::contains("mean")), ~ . > expthres))

    return(dfstrand)
}, unique(exptab$strand), unique(exptab$direction), MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)


exptranslist <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>% dplyr::arrange(transcript) %>% dplyr::pull(transcript) # nolint

return(list(main_table=alldf,exptranlist=exptranlist))
!!!!!!!!!!!!!!!


if (isTRUE(all.equal(length(idxscores), 0)))
    stop("The scores were not retrieved in the data.frame") # nolint
expressed_transcript_name <- main_table %>%
    group_by(transcript) %>%
    dplyr::summarize(gene=gene[1],strand=strand[1],
                    across(all_of(score_columns), ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint


# idxremove <- unique(unlist(apply(alldf[, idxscores], 2,
#     function(currentcol, thres) { return(which(currentcol < thres)) },
#     expthres, simplify = FALSE)))
