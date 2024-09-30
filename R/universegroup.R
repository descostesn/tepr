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
