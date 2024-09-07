########################
## This script aims at reproducing the figure 2d of the article using the
## object completedf generated with the script downstream.R
##
## Descostes - R-4.4.1 - sept 2024
########################

library("ggplot2")
options(digits = 2)


##################
# PARAMETERS
##################


completedfpath <- "/g/romebioinfo/tmp/downstream/completedf.rds"
dfmeandiffpath <- "/g/romebioinfo/tmp/downstream/dfmeandiff.rds"
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
colvec <- c("#90AFBB", "#FF9A04", "#10AFBB", "#FC4E07")
genename <- "EGFR"


##################
# MAIN
##################

## Reading inputs
completedf <- readRDS(completedfpath)
dfmeandiff <- readRDS(dfmeandiffpath)
expdf <- read.csv(exptabpath, header = TRUE)

plotecdf <- function(dfmeandiff, completedf, genename, verbose = TRUE) {

    ## Retrieving rows concerning the gene of interest
    if (verbose) message("\t Retrieving rows concerning the gene of interest")
    df <- dfmeandiff[which(dfmeandiff$gene == genename), ]

    ## Retrieving gene info
    geneinfo <- completedf[which(completedf$gene == genename), ]

    ## Computing the window size factor
    df100 <- df[which(df$coord == 100), ]
    windsizefact <- (df100$end - df100$start) / 1000


    !!!!!!!!!!

    AUC_cond_name <- paste0("AUC_", cond)
    AUC_cond <- round(gene_stat %>% pull(AUC_cond_name), decimal_places)  # Round AUC value
    KS_cond_name <-  paste0("adjFDR_p_AUC_", cond)
    KS_cond <- formatC(gene_stat %>% pull(KS_cond_name), format = "e", digits = decimal_places)
    knee_cond_name <- paste0("knee_AUC_", cond)
    knee_cond <- round(gene_stat %>% mutate(knee_kb=!!sym(knee_cond_name)*window_size/1000) %>% pull(knee_kb))  # Round Knee value
    

    # If KS test value is less than 0.01, store the knee_cond value for the vertical line
    
    if (as.numeric(KS_cond) < 0.01) {
      vline_data <- rbind(vline_data, data.frame(Condition = cond, Knee_Value = knee_cond))
      
    } else {
      knee_cond <- NA
    }
    
    sapply(unique(expdf$condition), function(cond) {
        return(paste0(geneinfo[,paste0("auc_", cond)], geneinfo[,paste0("adjFDR_", cond)]))
    })
    conditiontext <- paste0(cond, ": AUC = ", AUC_cond, ", KS = ", KS_cond, ", Knee (kb) = ", knee_cond)
    subtitle_text <- paste(subtitle_text, condition_text, sep = "\n")

    !!!!!!!!!!!!
}
