library("dplyr")
library("parallel")
library("tidyr")
library("tidyselect")

## /g/romebioinfo/tmp/comparewithscratch-downstream



##################
# PARAMETERS
##################


vicbigtsvpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/dTAG_Cugusi_stranded_20230810.tsv" # nolint
expdfpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
filtertabpath <- "/g/romebioinfo/Projects/tepr/Dataset/filtertab.csv"
nicbigtsvpath <- "/g/romebioinfo/tmp/comparewithscratch/finaltab.rds"
oldnictsvpath <- "/g/romebioinfo/tmp/preprocessing/backup/completeframedf.rds"
expthres <- 0.1
outputfolder <- "/g/romebioinfo/tmp/comparewithscratch"
nbcputrans <- 20
windsize <- 200



##################
#FUNCTIONS - downstream
##################














.retrievekneeandmax <- function(condvec, transtable) { # nolint

reslist <- lapply(condvec, function(cond, transtable) {
    difffxname <- paste0("diff_Fx_", cond)
    difffxvec <- transtable[, difffxname]
    ## If equality of difference within the same gene it takes the closest
    ## knee from the TSS # nolint
    resrow <- transtable[which(difffxvec == max(difffxvec)), ] %>% # nolint
        dplyr::slice_min(coord, n = 1) # nolint
    res <- data.frame(resrow$coord, resrow[, difffxname])
    colnames(res) <- c(paste0("knee_AUC_", cond), paste0("max_", difffxname))
    return(res)
}, transtable)

return(reslist)
}

kneeid <- function(transdflist, expdf, nbcputrans, verbose = FALSE) {

  condvec <- unique(expdf$condition)
  bytransres <- parallel::mclapply(transdflist, function(transtable, condvec) {
      bycondreslist <- .retrievekneeandmax(condvec, transtable)
       return(cbind(transcript = transtable$transcript[1],
        do.call("cbind", bycondreslist)))
    }, condvec, mc.cores = nbcputrans)
  res <- do.call("rbind", bytransres)
  return(res)
}


.summarytrans <- function(bytransmeanlist, nbcpu) {
  summarydflist <- mclapply(bytranslistmean, function(trans) {
    coor1 <- min(trans$coor1)
    coor2 <- max(trans$coor2)
    return(data.frame(chr = trans$chr[1], coor1, coor2,
          strand = trans$strand[1], gene = trans$gene[1],
          transcript = trans$transcript[1], size = coor2 - coor1 + 1))
  }, mc.cores = nbcpu)
  summarydf <- do.call("rbind", summarydflist)
  return(summarydf)
}

.computeupdown <- function(completbytrans, condvec, nbcpu) {

  updownbytranslist <- mclapply(completbytrans, function(trans, condvec) {

    ## Ordering by coordinates (security)
    trans <- trans[order(trans$coord), ]

    ## For each condition
    updownlist <- lapply(condvec, function(cond, trans) {
      kneecolname <- paste0("knee_AUC_", cond)
      meancolname <- paste0("mean_value_", cond)

      idxup <- which(trans$coord <= trans[, kneecolname])
      if (isTRUE(all.equal(length(idxup), 0)))
        stop("Problem in retrieving idxup, contact the developer.")
      upmean <- mean(trans[idxup, meancolname])

      idxdown <- which(trans$coord >= trans[, kneecolname] &
                          trans$coord <= max(trans$coord))
      if (isTRUE(all.equal(length(idxdown), 0)))
        stop("Problem in retrieving idxdown, contact the developer.")
      downmean <- mean(trans[idxdown, meancolname])

      ## Calculating attenuation
      att <- 100 - downmean / upmean * 100

      res <- data.frame(trans$transcript[1], upmean, downmean, att)
      colnames(res) <- c("transcript", paste0("UP_mean_", cond),
              paste0("DOWN_mean_", cond), paste0("Attenuation_", cond))
      return(res)
    }, trans)

    return(do.call("cbind", updownlist))
  }, condvec, mc.cores = nbcpu)

  updowndf <- do.call("rbind", updownbytranslist)
  updowndf <- updowndf[, -which(duplicated(colnames(updowndf)))]
  return(updowndf)
}

.filterattenuation <- function(auckneenasumatt, condvec, pval, replaceval, # nolint
    verbose) {

        mat <- auckneenasumatt
        if (verbose) message("\t\t Replacing non-significant attenuations by ",
            replaceval)
        invisible(sapply(condvec, function(cond, replaceval) {
            pauccond <- paste0("p_AUC_", cond)
            ## Replacing Attenuation value if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("Attenuation_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("Attenuation_", cond)]]))
            ## Replacing knee values if KS test > pval
            mat <<- mat %>%
                dplyr::mutate(!!paste0("knee_AUC_", cond) := # nolint
                    ifelse(.data[[pauccond]] >= pval, replaceval, # nolint
                    .data[[paste0("knee_AUC_", cond)]]))
        }, replaceval))

        return(mat)
}

attenuation <- function(allaucdf, kneedf, matnatrans, bytranslistmean, expdf,
  dfmeandiff, nbcpu = 1, significant = FALSE, replaceval = NA, pval = 0.1,
  verbose = TRUE) {

      if (verbose) message("\t Merging tables")
      allaucknee <- merge(allaucdf, kneedf, by = "transcript")
      mergecolnames <- c("gene", "transcript", "strand")
      allauckneena <- merge(allaucknee, matnatrans, by = mergecolnames)

      if (verbose) message("\t Building summary")
      summarydf <- .summarytrans(bytranslistmean, nbcpu)
      if (verbose) message("\t Merging summary")
      auckneenasum <- merge(summarydf, allauckneena, by = mergecolnames)

      ## Merging the mean table with the previous one
      if (verbose) message("\t Merging detailed mean table with summary")
      complet <- merge(dfmeandiff, auckneenasum, by = mergecolnames)

      ## Splitting the previous table by transcript
      if (verbose) message("\t Splitting the previous table by transcript")
      completbytrans <- split(complet, factor(complet$transcript))
      condvec <- unique(expdf$condition)

      ## For each transcript
      if (verbose) message("\t Computing up and down mean")
      updowndf <- .computeupdown(completbytrans, condvec, nbcpu)

      ## Merging attenuation to the complete table
      if (verbose) message("\t Merging attenuation to the complete table")
      auckneenasumatt <- merge(auckneenasum, updowndf, by = "transcript")

      ## Replace the attenuation values by replaceval if p_AUC_cond >= pval
      if (significant) {
        if (verbose) message("\t Keeping significant attenuation")
        auckneenasumatt <- .filterattenuation(auckneenasumatt, condvec, pval,
            replaceval, verbose)
      }

      return(auckneenasumatt)
}


universegroup <- function(completedf, controlname = "ctrl", stressname = "HS", # nolint
    windsizethres = 50, countnathres = 20, meanctrlthres = 0.5,
    meanstressthres = 0.5, pvaltheorythres = 0.1, aucctrlthreshigher = -10,
    aucctrlthreslower = 15, aucstressthres = 15, attenuatedpvalksthres = 2,
    outgrouppvalksthres = 0.2) {

    meanctrl <- paste("MeanValueFull", controlname, sep = "_")
    meanstress <- paste("MeanValueFull", stressname, sep = "_")
    pvaltheory <- paste("adjFDR_p_AUC", controlname, sep = "_")
    aucctrl <- paste("AUC", controlname, sep = "_")
    aucstress <- paste("AUC", stressname, sep = "_")
    pvalks <- paste0("adjFDR_p_dAUC_Diff_meanFx_", stressname, "_", controlname)

    ## Computing the Universe column
    completedf <- completedf %>%
        dplyr::mutate(Universe = ifelse(window_size > windsizethres & # nolint
            Count_NA < countnathres & !!sym(meanctrl) > meanctrlthres & # nolint
            !!sym(meanstress) > meanstressthres &
            !!sym(pvaltheory) > pvaltheorythres, TRUE, FALSE)) %>%
            dplyr::relocate(Universe, .before = 1)  # nolint

    ## Computing the Group column
    completedf <- completedf %>%
        dplyr::mutate(
            Group = ifelse(Universe == TRUE &
                !!sym(aucstress) > aucstressthres &
                -log10(!!sym(pvalks)) > attenuatedpvalksthres, "Attenuated",
                NA),
            Group = ifelse(Universe == TRUE &
                !!sym(pvalks) > outgrouppvalksthres &
                !!sym(aucctrl) > aucctrlthreshigher &
                !!sym(aucctrl) < aucctrlthreslower, "Outgroup", Group)) %>% # nolint
                dplyr::relocate(Group, .before = 2)

    return(completedf)
}



##################
# MAIN
##################


## This code tests the functions of downstream.R with the input table of
## victor: /g/romebioinfo/Projects/tepr/testfromscratch/bedgraph255/dTAG_Cugusi_stranded_20230810.tsv # nolint
bigtsv <- read.table(vicbigtsvpath, header = FALSE)
expdf <- read.csv(expdfpath, header = TRUE)


####
#### averageandfilterexprs
####

niccode_allexprsdfsvic <- averageandfilterexprs(expdf, bigtsv, expthres,
    verbose = TRUE)
viccode_allexprsdfsvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/results_main_table.rds") # nolint

if (isTRUE(all.equal(niccode_allexprsdfsvic[[1]], viccode_allexprsdfsvic[[1]])))
    message("table is equal after averageandfilterexprs")

if (isTRUE(all.equal(niccode_allexprsdfsvic[[2]], viccode_allexprsdfsvic[[2]])))
    message("transcript list is equal after averageandfilterexprs")

saveRDS(niccode_allexprsdfsvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_allexprsdfsvic.rds") # nolint

####
#### genesECDF
####

niccode_resecdfvic <- genesECDF(niccode_allexprsdfsvic, expdf, nbcpu = nbcputrans, verbose = TRUE) # nolint
nbwindows <- niccode_resecdfvic[[2]]
niccode_resecdfvic <- niccode_resecdfvic[[1]]

saveRDS(niccode_resecdfvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_resecdfvic.rds") # nolint

## Reading the result of ecdf that contains the column coord that is present in
## the input table of nic
viccode_resecdfvicpath <- "/g/romebioinfo/Projects/tepr/testfromscratch/cugusi2023_ECDFScores_10_200.tsv" # nolint
viccode_resecdfvic <- read.table(viccode_resecdfvicpath, sep = "\t", header = TRUE) # nolint

if (isTRUE(all.equal(as.data.frame(niccode_resecdfvic), viccode_resecdfvic)))
    message("genesECDF is consistent")


####
#### meandifference
####

## IMPORTANT: For the sake of comparison with the code of vic, only the first
## part of meandifference was executed by adding the following lines after
## resmean:
##        res <- cbind(resultsecdf, resmean)
##        return(res)

viccode_dfmeanvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/concat_dfFX_res.rds") # nolint
niccode_dfmeanvic <- meandifference(niccode_resecdfvic, expdf, nbwindows)

if (isTRUE(all.equal(viccode_dfmeanvic, niccode_dfmeanvic)))
    message("consistancy after dfmean")

## IMPORTANT: Now the whole function is executed (above lines are commented) to
## compute the differences of means
viccode_dfmeandiffvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/concat_Diff_mean_res.rds") # nolint
niccode_dfmeandiffvic <- meandifference(niccode_resecdfvic, expdf, nbwindows)

saveRDS(niccode_dfmeandiffvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_dfmeandiffvic.rds") # nolint

## Change the 'V' of 'value' to lower case in vic table
colnames(viccode_dfmeandiffvic)[which(colnames(viccode_dfmeandiffvic) == "Diff_meanValue_ctrl_HS")] <- "Diff_meanvalue_ctrl_HS" # nolint
colnames(viccode_dfmeandiffvic)[which(colnames(viccode_dfmeandiffvic) == "Diff_meanValue_HS_ctrl")] <- "Diff_meanvalue_HS_ctrl" # nolint

if (isTRUE(all.equal(viccode_dfmeandiffvic, niccode_dfmeandiffvic)))
    message("consistancy after mean differences")


####
#### dAUC
####

bytranslistmean <- split(niccode_dfmeandiffvic,
    factor(niccode_dfmeandiffvic$transcript))

## IMPORTANT: Only the delta auc (first part) of the function allauc is executed
## by adding return(daucallcond) in the function. This should be removed later.
niccode_daucdfvic <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)
viccode_daucdfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/dAUC_allcondi_res.rds") # nolint

# For comparison only
rownames(niccode_daucdfvic) <- NULL
names(viccode_daucdfvic$D_dAUC_Diff_meanFx_HS_ctrl) <- NULL

if (isTRUE(all.equal(viccode_daucdfvic, niccode_daucdfvic)))
    message("consistancy after dAUC")

####
#### AUC
####

## IMPORTANT: Only the second part of allauc is executed here by skipping the
## code when pasting in the terminal and adding return(aucallcond)
niccode_aucdfvic <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)
viccode_aucdfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/AUC_allcondi_res.rds") # nolint

# For comparison only
rownames(niccode_aucdfvic) <- rownames(viccode_aucdfvic) <- NULL
names(viccode_aucdfvic$D_AUC_ctrl) <- names(viccode_aucdfvic$D_AUC_HS) <- NULL

if (isTRUE(all.equal(viccode_aucdfvic, niccode_aucdfvic)))
    message("consistancy after AUC")


####
#### countNA
####

niccode_countnavic <- countna(niccode_allexprsdfsvic, expdf, nbcputrans)
viccode_countnavic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/count_NA_res.rds") # nolint

## For comparison only
niccode_countnavictest <- niccode_countnavic[order(niccode_countnavic$gene), ]
rownames(niccode_countnavictest) <- NULL
viccode_countnavictest <- as.data.frame(viccode_countnavic)
viccode_countnavictest <- viccode_countnavictest[order(viccode_countnavictest$gene), ] # nolint
rownames(viccode_countnavictest) <- NULL

if (isTRUE(all.equal(viccode_countnavictest, niccode_countnavictest)))
    message("consistancy after countna")


####
#### kneeid
####

bytranslistmean <- split(niccode_dfmeandiffvic,
    factor(niccode_dfmeandiffvic$transcript))

niccode_kneedfvic <- kneeid(bytranslistmean, expdf, nbcputrans)
viccode_kneedfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/KneeID_res.rds") # nolint

saveRDS(bytranslistmean, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/bytranslistmean.rds") # nolint
saveRDS(niccode_kneedfvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_kneedfvic.rds") # nolint

## for comparison only
rownames(niccode_kneedfvic) <- NULL

print(all.equal(niccode_kneedfvic, viccode_kneedfvic))

idxdiffctrl <- which(niccode_kneedfvic$knee_AUC_ctrl != viccode_kneedfvic$knee_AUC_ctrl) # nolint
idxdiffHS <- which(niccode_kneedfvic$knee_AUC_HS != viccode_kneedfvic$knee_AUC_HS) # nolint

message("The number of different knee in ctrl is: ", length(idxdiffctrl),
    ". The value for nic is ", niccode_kneedfvic[idxdiffctrl, "knee_AUC_ctrl"],
    " and the value for vic is ",
    viccode_kneedfvic[idxdiffctrl, "knee_AUC_ctrl"])
# The number of different knee in ctrl is: 1. The value for nic is 95 and the value for vic is 94 # nolint

message("The number of different knee in ctrl is: ", length(idxdiffHS),
    ". The value for nic is ", niccode_kneedfvic[idxdiffHS, "knee_AUC_HS"],
    " and the value for vic is ",
    viccode_kneedfvic[idxdiffHS, "knee_AUC_HS"])
# The number of different knee in ctrl is: 1. The value for nic is 65 and the value for vic is 64 # nolint



####
#### Attenuation
####


## Recomputing both types of auc
niccode_allaucdfvic <- allauc(bytranslistmean, expdf, nbwindows, nbcputrans)

saveRDS(niccode_allaucdfvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_allaucdfvic.rds") # nolint

niccode_completedfvic <- attenuation(niccode_allaucdfvic, niccode_kneedfvic,
    niccode_countnavic, bytranslistmean, expdf, niccode_dfmeandiffvic,
    nbcpu = nbcputrans)
backupniccode_completedfvic <- niccode_completedfvic

saveRDS(niccode_completedfvic, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_completedfvic.rds") # nolint

viccode_completedfvic <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/tst_df.rds") # nolint
viccode_completedfvic <- as.data.frame(viccode_completedfvic)

idx <- match(colnames(viccode_completedfvic), colnames(niccode_completedfvic))
niccode_completedfvic <- niccode_completedfvic[, idx]

## As highlighted previously, the knee values of rows 11512 and 6423 are
## different by a rank of 1. Therefore, these two rows are made equal for
## further comparison
niccode_completedfvic[c(11512, 6423), ] <- viccode_completedfvic[c(11512, 6423), ] # nolint

if (isTRUE(all.equal(niccode_completedfvic, viccode_completedfvic)))
    message("consistancy after attenuation")

## Testing the attenuation filtering
niccode_completedfvicfilt <- attenuation(niccode_allaucdfvic, niccode_kneedfvic,
    niccode_countnavic, bytranslistmean, expdf, niccode_dfmeandiffvic,
    nbcpu = nbcputrans, significant = TRUE)

viccode_completedfvicfilt <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/tst_dffilt.rds") # nolint
viccode_completedfvicfilt <- as.data.frame(viccode_completedfvicfilt)

idxfilt <- match(colnames(viccode_completedfvicfilt),
    colnames(niccode_completedfvicfilt))
niccode_completedfvicfilt <- niccode_completedfvicfilt[, idxfilt]

niccode_completedfvicfilt[c(11512, 6423), ] <- viccode_completedfvicfilt[c(11512, 6423), ] # nolint

if (isTRUE(all.equal(niccode_completedfvicfilt, viccode_completedfvicfilt)))
    message("consistancy after attenuation and filtering")

saveRDS(niccode_completedfvicfilt, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_completedfvicfilt.rds") # nolint

####
#### Universe and group
####

niccode_unigroupdf <- universegroup(niccode_completedfvic)
viccode_unigroupdf <- readRDS("/g/romebioinfo/Projects/tepr/testfromscratch/universegroupdf.rds") # nolint

if (isTRUE(all.equal(table(viccode_unigroupdf$Universe),
    table(niccode_unigroupdf$Universe))))
    message("Universe column is consistant")

if (isTRUE(all.equal(table(viccode_unigroupdf$Group),
    table(niccode_unigroupdf$Group))))
    message("Group column is consistant")

saveRDS(niccode_unigroupdf, file = "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_unigroupdf.rds") # nolint
