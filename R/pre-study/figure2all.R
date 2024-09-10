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
victabpath <- "/g/romebioinfo/Projects/tepr/downloads/Cugusi2022_AttenuationScores_10_200.tsv" # nolint
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
outputfolder <- "/g/romebioinfo/tmp/figures"


##################
#FUNCTIONS
##################

plotauc <- function(df, expdf, genevec, outfile = "AUCcompare_pval.pdf", # nolint
    outfold = "./", formatname = "pdf", plot = FALSE) {

    ## Retrieving the column containing the adj pval for delta acu ks.test
    colnamevec <- colnames(df)
    kstestadjname <- colnamevec[grep("adjFDR_pvaldeltadaucks_mean", colnamevec)]
    ## Retrieving the values in -log10 in a new column
    df <- cbind(df, kstestlog10 = -log10(df[, kstestadjname]))

    ## Building the names of the auc columns according to the conditions
    condvec <- unique(expdf$condition)
    idxctrl <- grep("ctrl", condvec)
    auccols <- paste("auc", condvec, sep = "_")
    ksvals <- "kstestlog10"

    ## Building data.frame of the genes to highlight
    genedf <- subset(df, gene %in% genevec) # nolint

    ## Structure of the basic scatterplot
    g <- ggplot2::ggplot(df %>% dplyr::arrange(df[, ksvals]),
        ggplot2::aes_string(auccols[idxctrl], auccols[-idxctrl],
        color = ksvals)) +
        ggplot2::geom_point(size = 0.5) + ggplot2::geom_density_2d()

    ## Adding highlight of the genes
    g1 <- g + ggrepel::geom_label_repel(data = genedf, aes(label = gene),
                                        box.padding   = 0.55,
                                        point.padding = 0,
                                        segment.color = 'black',
                                        max.overlaps = 50,
                                        color = "red")

    ## Formatting functions
    g2 <- g1 + ggplot2::scale_color_gradient2(midpoint = 0,  low = "white",
        mid = "grey", high = "darkgreen") +
        ggplot2::xlim(-10, 100) + ggplot2::ylim(-10, 100) +
        ggplot2::labs(x = "AUC in Control",
            y = paste0("AUC in ", condvec[-idxctrl]),
            legend = "-log10 p-value", color = "-log10 p-value") +
        ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "bottom")

    if (plot)
        print(g2)
    else
        ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
            plot = g2, device = formatname, path = outfold)
}


##################
# MAIN
##################

unigroupdf <- readRDS(unigrouppath)
expdf <- read.csv(exptabpath, header = TRUE)
genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")
victab <- read.delim(victabpath, header = TRUE)

## Plotting scatter of auc per condition
genevec <- c("EGFR","DAP","FLI1","MARCHF6", "LINC01619")
plotauc(unigroupdf, expdf, genevec, outfold = outputfolder, plot = TRUE)

## Test plot on vic tab

plotaucvic <- function()
auc_ctrlname <- "AUC_ctrl"
auc_stressname <- "AUC_HS"
pvalkstest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
df <- cbind(victab, kstestlog10 = -log10(victab[, pvalkstest]))
kstestlog10str <- "kstestlog10"
labelx <- "AUC in Control"
labely <- "AUC in HS"
axismin_x <- -10
axismax_x <- 100
axismin_y <- -10
axismax_y <- 100

g <- ggplot2::ggplot(df %>% dplyr::arrange(df[, kstestlog10str]),
    ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname), color= !!sym(kstestlog10str))) +
  ggplot2::geom_point(size=0.5) + ggplot2::geom_density_2d()

g1 <- g + ggrepel::geom_label_repel(data = subset(df, gene %in% genevec), aes(label = gene),
   box.padding   = 0.55,
   point.padding = 0,
   segment.color = 'black', max.overlaps = 50, color="red")

g2 <- g1 +  ggplot2::scale_color_gradient2(midpoint=0,  low="white", mid="grey", high = "darkgreen") +
  ggplot2::xlim(axismin_x,axismax_x) + ggplot2::ylim(axismin_y,axismax_y)+
  ggplot2::labs(x=labelx, y=labely, legend="-log10 p-value", color="-log10 p-value") +
  ggplot2::coord_fixed(ratio = 1) +   # Set aspect ratio to 1:1
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "bottom" )


!!!!!!!!!

auccols <- c("AUC_ctrl", "AUC_HS")
idxctrl <- 1
condvec <- unique(expdf$condition)
ksvals <- "p_dAUC_Diff_meanFx_HS_ctrl"
genedf <- subset(victab, gene %in% genevec) # nolint
g <- ggplot2::ggplot(victab %>% dplyr::arrange(victab[, ksvals]),
        ggplot2::aes_string(auccols[idxctrl], auccols[-idxctrl],
        color = ksvals)) +
        ggplot2::geom_point(size = 0.5) + ggplot2::geom_density_2d()

    ## Adding highlight of the genes
    g1 <- g + ggrepel::geom_label_repel(data = genedf, aes(label = gene),
                                        box.padding   = 0.55,
                                        point.padding = 0,
                                        segment.color = 'black',
                                        max.overlaps = 50,
                                        color = "red")

    ## Formatting functions
    g2 <- g1 + ggplot2::scale_color_gradient2(midpoint = 0,  low = "white",
        mid = "grey", high = "darkgreen") +
        ggplot2::xlim(-10, 100) + ggplot2::ylim(-10, 100) +
        ggplot2::labs(x = "AUC in Control",
            y = paste0("AUC in ", condvec[-idxctrl]),
            legend = "-log10 p-value", color = "-log10 p-value") +
        ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "bottom")
