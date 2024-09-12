########################
## This script aims at reproducing the panels of figure 2 (excluding panel d)
## of the article using the objects XX generated with the script downstream.R
##
## Descostes - R-4.4.1 - sept 2024
########################


library("ggplot2")
library("dplyr")
library("ggrepel")


##################
# PARAMETERS
##################

unigrouppath <- "/g/romebioinfo/tmp/downstream/unigroupdf.rds"
exploretabpath <- "/g/romebioinfo/tmp/explore/tst_df-3.rds"
victabpath <- "/g/romebioinfo/Projects/tepr/downloads/Cugusi2022_AttenuationScores_10_200.tsv" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
outputfolder <- "/g/romebioinfo/tmp/figures"
genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")



##################
#FUNCTIONS
##################

plotauc <- function(tab, auc_ctrlname, auc_stressname, pvalkstestcolname, # nolint
    labelx = "AUC in Control", labely = "AUC in Stress", axismin_x = -10,
    axismax_x = 100, axismin_y = -10, axismax_y = 100, maintitle = "",
    subtitle = "", legendpos = "bottom", formatname = "pdf", outfold = "./",
    outfile = "AUCcompare_pval.pdf", plottype = "pval", plot = FALSE,
    universename = "universe", groupname = "group") {

        if (!isTRUE(all.equal(plottype, "pval")) &&
            !isTRUE(all.equal(plottype, "groups")))
            stop("plottype should be equal to 'pval' or 'groups'.")

        if (isTRUE(all.equal(plottype, "pval"))) {
            df <- cbind(tab, kstestlog10 = -log10(tab[, pvalkstestcolname]))
            kstestlog10str <- "kstestlog10"
            df <- df %>% dplyr::arrange(df[, kstestlog10str])
            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname), # nolint
             color = !!sym(kstestlog10str))
            geompointinfo <- ggplot2::geom_point(size = 0.5)
            geompointinfo2 <- ggplot2::geom_density_2d()
        } else {
            df <- tab %>% dplyr::filter(!!sym(universename)) # nolint
            dfatt <- tab %>% dplyr::filter(!!sym(groupname) == "Attenuated") # nolint
            dfoutgroup <- tab %>% filter(!!sym(groupname) == "Outgroup") # nolint

            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname)) # nolint
            geompointinfo <- ggplot2::geom_point(size = 0.5, color = "grey")
            geompointinfo2 <- ggplot2::geom_point(data = dfatt, aesvar,
                color = "#e76f51", size = 1)
            geompointinfo3 <- ggplot2::geom_point(data = dfoutgroup, aesvar,
                    color = "#e9c46a", size = 1)
        }

        ## Structure of the basic scatterplot
        g <- ggplot2::ggplot(df, aesvar) + geompointinfo + geompointinfo2

        if (isTRUE(all.equal(plottype, "pval"))) {
            ## Adding highlight of the genes
            g <- g + ggrepel::geom_label_repel(data = subset(df,
                gene %in% genevec), aes(label = gene), box.padding = 0.55, # nolint
                point.padding = 0, segment.color = "black", max.overlaps = 50,
                color = "red") +
                ggplot2::scale_color_gradient2(midpoint = 0, low = "white",
                    mid = "grey", high = "darkgreen")
        } else {
            g <- g + geompointinfo3
        }

        ## Formatting functions
        g <- g + ggplot2::xlim(axismin_x, axismax_x) +
            ggplot2::ylim(axismin_y, axismax_y) +
            ggplot2::labs(x = labelx, y = labely, legend = "-log10 p-value",
                color = "-log10 p-value", title = maintitle,
                subtitle = subtitle) +
            ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_classic() +
            ggplot2::theme(legend.position = legendpos)

        if (plot) {
            warning("You chose to plot the auc, the figure is not saved.")
            print(g)
        } else {
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}

victorpreprocess <- function(victab) {

    mean_value_control_full <- "MeanValueFull_ctrl"
    mean_value_stress <- "MeanValueFull_HS"
    AUC_ctrl <- "AUC_ctrl"
    AUC_stress <- "AUC_HS"
    p_value_KStest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
    p_value_theoritical<- "adjFDR_p_AUC_ctrl" 

    tst_df <- victab %>%
    mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 &
        !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
        !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
        relocate(Universe, .before = 1) 
    tst_df <- tst_df %>%
    mutate(Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 &
    -log10(!!sym(p_value_KStest))>2, "Attenuated", NA),
        Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 &
        !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group)
    ) %>%  relocate(Group, .before = 2)

    return(tst_df)
}



##################
# MAIN
##################

unigroupdf <- readRDS(unigrouppath)
expdf <- read.csv(exptabpath, header = TRUE)
victab <- read.delim(victabpath, header = TRUE)
exploretab <- readRDS("/g/romebioinfo/tmp/explore/tst_df-3.rds")

## Performing preprocessing on victab - REMOVE FROM FINAL CODE
tst_df <- victorpreprocess(victab)

## Test plot with pval on vic tab - REMOVE FROM FINAL CODE
plotauc(tst_df, "AUC_ctrl", "AUC_HS", "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in HS", outfold = outputfolder,
    plot = TRUE)

## Test plot with pval on nic tab
plotauc(unigroupdf, "auc_ctrl", "auc_HS", "adjFDR_pvaldeltadaucks_mean_Fx_HS",
    labelx = "AUC in Control", labely = "AUC in HS", outfold = outputfolder,
    plot = TRUE)

## Test plot with groups on vic tab - REMOVE FROM FINAL CODE
plotauc(tst_df, "AUC_ctrl", "AUC_HS", "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in HS", legendpos = "none",
    subtitle = "Genes selected for Unibind", maintitle = "AUC Control vs HS",
    outfold = outputfolder, outfile = "AUCcompare_groups.pdf", plot = TRUE,
    plottype = "groups", universename = "Universe", groupname = "Group")
plotauc(exploretab, "AUC_ctrl", "AUC_HS", "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in HS", legendpos = "none",
    subtitle = "Genes selected for Unibind", maintitle = "AUC Control vs HS",
    outfold = outputfolder, outfile = "AUCcompare_groups.pdf", plot = TRUE,
    plottype = "groups", universename = "Universe", groupname = "Group")

## Test plot with groups on nic tab
plotauc(unigroupdf, "auc_ctrl", "auc_HS", "adjFDR_pvaldeltadaucks_mean_Fx_HS",
    labelx = "AUC in Control", labely = "AUC in HS", legendpos = "none",
    subtitle = "Genes selected for Unibind", maintitle = "AUC Control vs HS",
    outfold = outputfolder, outfile = "AUCcompare_groups.pdf", plot = TRUE,
    plottype = "groups")
