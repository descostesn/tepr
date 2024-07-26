####################################
# This script goes through documentation/explore.R and homogenizes it with
# preprocessing.R
#
# Descostes - R-4.4.1 - July 2024
####################################

library("tidyr")
library("dplyr")
library("tidyselect")
library("parallel")


##################
# PARAMETERS
##################

alldfpath <- "/g/romebioinfo/Projects/tepr/robjsave/alldffrompreprocessing.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv"
expthres <- 0.1
nbcpu <- 10


##################
#FUNCTIONS
##################

averageandfilterexprs <- function(expdf, alldf, expthres) { # nolint

    scorecolvec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ## Calculate the average expression per transcript (over each frame)
    dfbytranscript <- alldf %>% dplyr::group_by(transcript) %>% # nolint
        dplyr::summarize(gene = gene[1], strand = strand[1], # nolint
            dplyr::across(
                tidyselect::all_of(scorecolvec),
                ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint
    ## Remove a line if it contains only values < expthres (separating strands)
    dfstrandlist <- mapply(function(strandname, directname, dfbytrans,
        expthres) {
            if ((isTRUE(all.equal(strandname, "+")) &&
                isTRUE(all.equal(directname, "rev"))) ||
                (isTRUE(all.equal(strandname, "-")) &&
                isTRUE(all.equal(directname, "fwd"))))
                    stop("Strand and direction do not match, contact the ",
                        "developer")
            dfstrand <- dfbytranscript %>%
                dplyr::filter(strand == strandname) %>% # nolint
                dplyr::select(gene, transcript, strand, # nolint
                tidyselect::contains(directname))  %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ !is.na(.))) %>%
                dplyr::filter(dplyr::if_all(tidyselect::all_of(
                tidyselect::contains("mean")), ~ . > expthres))
            return(dfstrand)
        }, unique(expdf$strand), unique(expdf$direction),
            MoreArgs = list(dfbytranscript, expthres), SIMPLIFY = FALSE)

    exptranstab <- dplyr::bind_rows(dfstrandlist[[1]], dfstrandlist[[2]]) %>%
        dplyr::arrange(transcript) %>% dplyr::pull(transcript) # nolint

    return(list(maintable = alldf, exptranlist = exptranstab))
}

.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname,
                " should be unique, contact the developer.")
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1) { # nolint

    ## Defining variables
    maintable <- allexprsdfs[[1]]
    exprstransnames <- allexprsdfs[[2]]

    ## Filtering the main table to keep only the expressed transcripts
    idx <- match(maintable$transcript, exprstransnames)
    idxnoexpr <- which(is.na(idx))
    maintable <- maintable[-idxnoexpr, ]

    ## Splitting the table by each transcript to perform transcript specific
    ## operations
    transdflist <- split(maintable, factor(maintable$transcript))
    nbrows <- unique(sapply(transdflist, nrow)) ## all transcripts have the same number of windows, no need to calculate it each time # nolint
    .checkunique(nbrows, "nbrows")
    framevec <- seq_len(nbrows)
    colnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        framevec, colnamevec) {
        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        .checkunique(str, "str")
        colnamestr <- colnamevec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- apply(scoremat, 2, function(x, rounding, framevec) {
            extendedframevec <- rep(framevec, ceiling(x * rounding))
            fx <- ecdf(extendedframevec)(framevec)
            return(fx)
        }, rounding, framevec, simplify = TRUE)
        colnames(ecdfmat) <- paste("Fx", colnames(ecdfmat), sep = "_")
        res <- cbind(transtable, ecdfmat)
        return(res)
    }, expdf, framevec, colnamevec, mc.cores = nbcpu)

    concatdf <- dplyr::bind_rows(ecdflist)

    return(concatdf)
}

meananddiff <- function(resultsecdf, exptab) {

    lapply(exptab$condition, function(currentcond, df) {

        message("Merging columns for condition ", currentcond)
        ## Retrieving columns having condition name as substring
        idxcond <- grep(currentcond, colnames(df))
        if (isTRUE(all.equal(length(idxcond), 0)))
            stop("Problem in function meananddiff, condition not found in ",
                "column names. Contact the developer.")

        ## Separating column names by values and Fx
        idxcondfx <- grep("Fx", colnames(df[idxcond]))
        if (isTRUE(all.equal(length(idxcondfx), 0)))
            stop("Problem in function meananddiff, column Fx not found in ",
                "column names. Contact the developer.")
        idxcondlist <- list(condval = idxcond[idxcondfx],
            condfx = idxcond[-idxcondfx])

        ## Computing mean and diff values for score and fx columns
    }, resultsecdf)
    
    
        # Calculate row means for the specified columns
        # Check if there is more than one replicate
        if (length(replicate_numbers) > 1) {
            ## !! I do not get how this can work since the column are not defined!!
          concat_df[[mean_value_condi_name]] <- rowMeans(concat_df[, column_vector_value], na.rm = F) # nolint
          concat_df[[mean_Fx_condi_name]] <- rowMeans(concat_df[, column_vector_Fx], na.rm = FALSE) # nolint
        } else { ## ! this case is not necessary, 
          # Handle case when column_vector_value is empty
          new_column_value <- paste0("value_", cond, "_rep", "1", "_score") # Generate a new item # nolint
          new_column_Fx <- paste0("Fx_", cond, "_rep", "1", "_score") # Generate a new item # nolint
## !! I do not get how this can work since the column are not defined!!
        concat_df[[mean_value_condi_name]] <- concat_df[[new_column_value]]
        concat_df[[mean_Fx_condi_name]] <- concat_df[[new_column_Fx]]
        }
## !! I do not get how this can work since the column are not defined!!
        concat_df[[diff_Fx_condi_name]] <- concat_df[[mean_Fx_condi_name]] - concat_df$coord/window_number ## Difference with the y=x ECDF, used to calculate AUC # nolint
        column_vector_value <- character() ## obligatory to reset the columns values to empty # nolint
        column_vector_Fx <- character() # nolint
    }

      return(concat_dfFx=concat_df)
}


##################
# MAIN
##################

## Reading alldf and info tab
alldf <- readRDS(alldfpath)
expdf <- read.csv(exptabpath, header = TRUE)

## Filtering out non expressed transcripts:
## 1) for each column, calculate the average expression per transcript (over each frame) # nolint
## 2) For each column, remove a line if it contains only values < expthres separating strands # nolint
allexprsdfs <- averageandfilterexprs(expdf, alldf, expthres)
resultsecdf <- genesECDF(allexprsdfs, expdf, nbcpu = nbcpu)
