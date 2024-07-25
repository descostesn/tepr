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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

.checkunique <- function(x, xname) {
        if (!isTRUE(all.equal(length(x), 1)))
            stop("The element ", xname,
                " should be unique, contact the developer.")
}

genesECDF <- function(allexprsdfs, expdf, rounding = 10, nbcpu = 1) {

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
    colnamevec <- paste0(expdf$condition, expdf$replicate, expdf$direction)

    ecdflist <- parallel::mclapply(transdflist, function(transtable, expdf,
        nbrows, colnamevec) {
        ## Filters the score columns according to the strand of the transcript
        str <- as.character(unique(transtable$strand))
        .checkunique(str, "str")
        colnamestr <- colnamevec[which(expdf$strand == str)]
        scoremat <- transtable[, colnamestr]

        ## For each column of the scoremat, compute ecdf
        ecdfmat <- apply(scoremat, 2, function(x, rounding, nbrows) {
            coordvec <- rep(seq_len(nbrows), ceiling(x*rounding))
            x <- x[coordvec]
            ecdfdf <- ecdf
        }, rounding, nbrows, simplify = TRUE)






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## This for loop is equivalent to computing on each column, bigDF might
        ## not be useful.
        for (my_var in unique(df_long$variable)) {
            df_subset <- subset(df_long, subset = variable == my_var) # This is just selecting the lines that we had in the initial table # nolint
            df_expanded <- df_subset[rep(seq_len(nrow(df_subset)), df_subset$value_round), ] # nolint
            ecdf_df <- ecdf(df_expanded[,"coord"])
            df_subset$Fx <- ecdf_df(df_subset$coord) 
            list_df[[i]] <- df_subset
            i <- i + 1
        }

        ## This two lines are equivalent to cbind a matrix of Fx to the transcript table used at the beginning
        ## The columns are biotype, chr, coor1, coor2, transcript, gene, strand, window, id, ctrl_rep1.plus, ctrl_rep2.plus,
        ## HS_rep1.plus, HS_rep2.plus, coord, value_ctrl_rep1.plus_score, value_ctrl_rep2.plus_score,
        ## value_HS_rep1.plus_score, value_HS_rep2.plus_score, Fx_ctrl_rep1.plus_score, Fx_ctrl_rep2.plus_score,
        ## Fx_HS_rep1.plus_score, Fx_HS_rep2.plus_score
        df_final <- bind_rows(list_df)
        transcript_table <- df_final  %>% pivot_wider(., names_from = "variable", values_from = c("value", "value_round", "Fx")) %>% select(., -contains("value_round")) # nolint
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }, expdf, nbrows, colnamevec, mc.cores = nbcpu)


    # res <- getting_var_names(extension, workdir)
    # col_names <- res$col_names
    # var_names <- res$var_names
    # total_iterations <- length(exprstransnames)
    # #setting the progress bar
    # pb <- txtProgressBar(min = 0, max = total_iterations, style =5)
    # j = 0
    # concat_df <- data.frame()
    ## Looping through all the transcripts i.e performs the loop for each
    ## transcript that has been kept by the function maintable_read because
    ## it was expressed
    # for (variable in exprstransnames) {
    #     gene_table <- data.frame()
    #     bigDF <- data.frame()

        ## Isolating the rows corresponding to the transcript
        ## Question: Could we do it at once for all expressed transcripts with
        ## a match (transcriptmaintable in selectedtrans, remove na, use the
        ## transcript column as factor to split the main table)
        # transcript_table <- data.frame()
        # transcript <- filter(maintable, maintable$transcript == variable)
        ## This replacement is not necessary in my hands
        # transcript[transcript == "NAN"] <- NA
        ## Changing the name of transcript can be useful to conserve it. So far,
        ## it seems to be done to keep looking for the strand. It could be
        ## stored in a variable. Add verification that the strand is unique
        ## in transcript.
        bigDF <- transcript
        ## This is equivalent to nrow(transcript). Check if every transcript has
        ## the same number of windows, if it is the case, it is not necessary
        ## to take this into account.
        my_length <- length(bigDF[,'window'])

# ---------------------------------------------------------------------------------------------------
## SUMMARY
##
## This section filters the score columns according to the strand of the transcript considered

        ## This only builds the names of the columns containing the scores
        var_names_score <- paste0(var_names,"_score")

        if (transcript$strand[1] == "-") { # see above, the strand can be stored in a variable instead of calculating it. # nolint
            bigDF <- bigDF %>%
            select(!matches("plus"))
            bigDF$coord <- seq(from = my_length, to = 1, by = -1)
            bigDF <- arrange(bigDF, coord)
            conditions <- var_names_score[grepl("minus", var_names_score)]
        } else {
            bigDF <- bigDF %>%
            select(!matches("minus"))
            bigDF$coord <- seq(from=1, to=my_length,by=1)
            conditions <- var_names_score[grepl("plus",var_names_score)]
        }
# ---------------------------------------------------------------------------------------------------

        ## To my understanding this Fills missing values in selected columns
        ## using first down and then up previous entry.
        bigDF <- bigDF %>% fill(contains("score"), .direction = "downup")

        ## The code below creates a data.frame of two columns with the name of
        ## the experiment as key and the expression as value. This is the kind
        ## of transformation we do before using ggplot2. However I do not see
        ## why keeping the previous columns is necessary. They will be repeated
        ## several times.
        df_long <- bigDF %>% 
            gather(key = "variable", value = "value", conditions)
        ## Makes the values as integer by multiplying by 10. Why not using the
        ## ceil or floor function. Is it important to increase the scores?
        df_long[,'value'] <- as.numeric(df_long[,'value'])
        df_long[,'value_round']<- round(df_long$value*rounding)

        #  Update the progress bar
        j = j + 1
        setTxtProgressBar(pb, j)

        list_df <- list()
        i <- 1
        ## This for loop is equivalent to computing on each column, bigDF might
        ## not be useful.
        for (my_var in unique(df_long$variable)) {
            df_subset <- subset(df_long, subset = variable == my_var) # This is just selecting the lines that we had in the initial table # nolint
            df_expanded <- df_subset[rep(seq_len(nrow(df_subset)), df_subset$value_round), ] # nolint
            ecdf_df <- ecdf(df_expanded[,"coord"])
            df_subset$Fx <- ecdf_df(df_subset$coord) 
            list_df[[i]] <- df_subset
            i <- i + 1
        }

        ## This two lines are equivalent to cbind a matrix of Fx to the transcript table used at the beginning
        ## The columns are biotype, chr, coor1, coor2, transcript, gene, strand, window, id, ctrl_rep1.plus, ctrl_rep2.plus,
        ## HS_rep1.plus, HS_rep2.plus, coord, value_ctrl_rep1.plus_score, value_ctrl_rep2.plus_score,
        ## value_HS_rep1.plus_score, value_HS_rep2.plus_score, Fx_ctrl_rep1.plus_score, Fx_ctrl_rep2.plus_score,
        ## Fx_HS_rep1.plus_score, Fx_HS_rep2.plus_score
        df_final <- bind_rows(list_df)
        transcript_table <- df_final  %>% pivot_wider(., names_from = "variable", values_from = c("value", "value_round", "Fx")) %>% select(., -contains("value_round")) # nolint

#---------------
## ADDED BY ME
modified_dataset <- transcript_table
#---------------


# ------------------------------------------------
## THIS CANNOT BE TRIGGERED BECAUSE THE COLUMNS DO NOT EXIST
        # getting rid of plus and minus
        if (transcript_table$strand[1]=="-") {
            # Drop columns containing "minus"
            columns_to_drop <- grep("plus", col_names, value = TRUE)
            dataset_without_dropped <- transcript_table %>%
            select(-all_of(columns_to_drop))

        # Modify column names by removing "_plus"
        modified_dataset <- dataset_without_dropped %>%
        rename_with(~gsub(".minus", "", .), contains(".minus"))
        } else {
            # Drop columns containing "minus"
            columns_to_drop <- grep("minus", col_names, value = TRUE)
            dataset_without_dropped <- transcript_table %>%
            select(-all_of(columns_to_drop))

            # Modify column names by removing "_plus"
            modified_dataset <- dataset_without_dropped %>%
            rename_with(~gsub(".plus", "", .), contains(".plus"))
        }
# ----------------------------------------------------
        concat_df <- bind_rows(concat_df, modified_dataset)
    }
saveRDS(concat_df, file = file.path("/g/romebioinfo/Projects/tepr/robjsave/concatdf_fromexplore.rds"))
    # # Close the progress bar
    close(pb)
    # list_gene_table <- concat_df %>% select(gene) %>% distinct()
    gc()

    return(concat_df = concat_df)
}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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




resultsECDF <- genesECDF(allexprsdfs, expdf, nbcpu = nbcpu)
