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
exptabpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab.csv" # nolint
outputfolder <- "/g/romebioinfo/tmp/figures"


##################
#FUNCTIONS
##################

plotauc <- function(df, expdf, genevec, outfile = "AUCcompare_pval.pdf", # nolint
    outfold = "./", formatname = "pdf") {

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
    g <- ggplot2::ggplot(df %>% dplyr::arrange(df[, "kstestlog10"]),
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

    ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
        plot = g2, device = formatname, path = outfold)
}


##################
# MAIN
##################

unigroupdf <- readRDS(unigrouppath)
expdf <- read.csv(exptabpath, header = TRUE)
genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")

## Plotting scatter of auc per condition
plotauc(unigroupdf, expdf, genevec, outfold = outputfolder)
