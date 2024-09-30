library("ggplot2")
library("tidyr")
library("tidyselect")

## /g/romebioinfo/tmp/comparewithscratch-plotting


##################
# PARAMETERS
##################

unigroupdfpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_unigroupdf.rds" # nolint
dfmeandiffpath <- "/g/romebioinfo/tmp/comparewithscratch-downstream/niccode_dfmeandiffvic.rds" # nolint
expdfpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/exptab-bedgraph-vicnames.csv" # nolint
colvec <- c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
outfold <- "/g/romebioinfo/tmp/comparewithscratch-plotting"


##################
#FUNCTIONS
##################




.callggplotauc <- function(df, aesvar, geompointinfo, geompointinfo2,
    geompointinfo3, plottype, axismin_x, axismax_x, axismin_y, axismax_y,
    labelx, labely, maintitle, subtitle, legendpos, plot, outfile, formatname,
    outfold, genevec) {

        ## Structure of the basic scatterplot
        g <- ggplot2::ggplot(df, aesvar) + geompointinfo + geompointinfo2

        if (isTRUE(all.equal(plottype, "pval")) && !is.na(genevec)) {

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
            ggplot2::labs(x = labelx, y = labely, legend = "-log10 p-value", # nolint
                color = "-log10 p-value", title = maintitle, # nolint
                subtitle = subtitle) +
            ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_classic() +
            ggplot2::theme(legend.position = legendpos)

        if (plot) {
            warning("You chose to plot the auc, the figure is not saved.") # nolint
            print(g)
        } else {
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}

.colnamecheck <- function(colnamevec, tab) {
            invisible(sapply(colnamevec, function(currentcol, tab) {
            idx <- grep(currentcol, colnames(tab))
            if (isTRUE(all.equal(length(idx), 0)))
                stop("The column ", currentcol, " does not exist in the ",
                    "provided table.")
        }, tab))
}

.checkplotaucparams <- function(plottype, auc_ctrlname, auc_stressname,
    pvalkstestcolname, genevec, tab) {

        if (!isTRUE(all.equal(plottype, "pval")) &&
            !isTRUE(all.equal(plottype, "groups")))
                stop("plottype should be equal to 'pval' or 'groups'.")

        colnamevec <- c(auc_ctrlname, auc_stressname, pvalkstestcolname)
        .colnamecheck(colnamevec, tab)

        if (isTRUE(all.equal(plottype, "groups")) && !is.na(genevec[1]))
            stop("The vector of genes is not necessary for plotting groups")
}

plotauc <- function(tab, genevec = NA, # nolint
    auc_ctrlname = "AUC_ctrl", auc_stressname = "AUC_HS",
    pvalkstestcolname = "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl",
    labelx = "AUC in Control", labely = "AUC in Stress", axismin_x = -10,
    axismax_x = 100, axismin_y = -10, axismax_y = 100, maintitle = "",
    subtitle = "", legendpos = "bottom", formatname = "pdf", outfold = "./",
    outfile = "AUCcompare_pval.pdf", plottype = "pval", plot = FALSE,
    universename = "Universe", groupname = "Group") {

        .checkplotaucparams(plottype, auc_ctrlname, auc_stressname,
            pvalkstestcolname, genevec, tab)

        if (isTRUE(all.equal(plottype, "pval"))) {
            df <- cbind(tab, kstestlog10 = -log10(tab[, pvalkstestcolname]))
            kstestlog10str <- "kstestlog10"
            df <- df %>% dplyr::arrange(df[, kstestlog10str])
            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname), # nolint
             color = !!sym(kstestlog10str))
            geompointinfo <- ggplot2::geom_point(size = 0.5)
            geompointinfo2 <- ggplot2::geom_density_2d()
        } else {
            df <- tab %>% dplyr::filter(!!sym(universename) == FALSE) # nolint
            dfatt <- tab %>% dplyr::filter(!!sym(groupname) == "Attenuated") # nolint
            dfoutgroup <- tab %>% dplyr::filter(!!sym(groupname) == "Outgroup") # nolint

            aesvar <- ggplot2::aes(!!sym(auc_ctrlname), !!sym(auc_stressname)) # nolint
            geompointinfo <- ggplot2::geom_point(size = 0.5, color = "grey")
            geompointinfo2 <- ggplot2::geom_point(data = dfatt, aesvar,
                color = "#e76f51", size = 1)
            geompointinfo3 <- ggplot2::geom_point(data = dfoutgroup, aesvar,
                    color = "#e9c46a", size = 1)
        }

        .callggplotauc(df, aesvar, geompointinfo, geompointinfo2,
            geompointinfo3, plottype, axismin_x, axismax_x, axismin_y,
            axismax_y, labelx, labely, maintitle, subtitle, legendpos, plot,
            outfile, formatname, outfold)
}


.normalizeandsummarize <- function(transvec, dfmeandiff, unigroupdf, daucname, # nolint
    auc_ctrlname, auc_stressname) {

    ## Selecting full mean and AUC columns
    AUC_allcondi <- unigroupdf %>% dplyr::select(transcript, gene, strand,
        dplyr::contains("Full"), !!sym(daucname), !!sym(auc_ctrlname),
        !!sym(auc_stressname), -contains(c("UP", "DOWN")), window_size)

    ## Selecting coord and mean values
    result <- dfmeandiff %>%
        dplyr::filter(transcript %in% transvec) %>%
        dplyr::left_join(., AUC_allcondi, by = c("transcript", "gene")) %>%
        dplyr::select(transcript, gene, coord, contains("mean_value"),
        -contains("Full"))  %>% dplyr::group_by(coord) %>%
        dplyr::summarise(dplyr::across(contains("mean_value"),
        ~ mean(., na.rm = TRUE)))

    return(result)
}

plotmetagenes <- function(unigroupdf, plottype = "attenuation",
    daucname = "dAUC_Diff_meanFx_HS_ctrl", auc_ctrlname = "AUC_ctrl",
    auc_stressname = "AUC_HS", plot = FALSE, formatname = "pdf",
    outfold = "./") {

    if (!isTRUE(all.equal(plottype, "attenuation")) &&
        !isTRUE(all.equal(plottype, "outgroup")) &&
        !isTRUE(all.equal(plottype, "universe")) &&
        !isTRUE(all.equal(plottype, "all")))
        stop("plot type should be one of: attenuation, outgroup, universe, all")

    colnamevec <- c(daucname, auc_ctrlname, auc_stressname)
    .colnamecheck(colnamevec, unigroupdf)

    ## Selection of transcripts and define plot title
    if (isTRUE(all.equal(plottype, "attenuation"))) {
        idx <- which(unigroupdf$Group == "Attenuated")
        titleplot <- "Attenuated genes"
    } else if (isTRUE(all.equal(plottype, "outgroup"))) {
        idx <- which(unigroupdf$Group == "Outgroup")
        titleplot <- "Outgroup genes"
    } else if (isTRUE(all.equal(plottype, "universe"))) {
        idx <- which(unigroupdf$Universe)
        titleplot <- "Universe genes"
    } else {
        idx <- seq_len(nrow(unigroupdf))
        titleplot <- "All genes"
    }

    if (isTRUE(all.equal(length(idx), 0)))
        stop("No transcripts were found for the criteria ", plottype)

     transvec <- unigroupdf[idx, "transcript"]
     df <- .normalizeandsummarize(transvec, dfmeandiff, unigroupdf, daucname,
        auc_ctrlname, auc_stressname)
    meanvalctrl <-  colnames(df)[2]
    meanvalstress <- colnames(df)[3]

    ## plotting
    g <-  ggplot2::ggplot() +
        ggplot2::geom_line(data = df, aes(x = coord/2,
        y = !!sym(meanvalctrl)), color = "#00AFBB", size = 1.5) +
        ggplot2::geom_line(data = df,
            aes(x = coord/2, y = !!sym(meanvalstress)),
            color = "#FC4E07", size = 1.5) +
        ggplot2::theme_bw() + ggplot2::ylim(0,7) +
        ggplot2::labs(x = "TSS to TTS", title = titleplot,
            subtitle = length(transvec), y = "Transcription density") +
        ggplot2::theme(legend.position = "none", legend.box = "vertical")

    if (plot) {
        warning("You chose to plot the auc, the figure is not saved.") # nolint
        print(g)
    } else {
        outfile <- paste0("metagene_", plottype)
        ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}


plothistoknee <- function(unigroupdf, plottype = "percent", xlimvec = NA, # nolint
    binwidthval = NA, kneename = "knee_AUC_HS", plot = FALSE, outfold = "./",
    formatname = "pdf", universename = "Universe", groupname = "Group") {

        if (!isTRUE(all.equal(plottype, "percent")) &&
            !isTRUE(all.equal(plottype, "kb")))
                stop("Plot type should be percent or kb")

        colnamevec <- c(universename, groupname, kneename)
        .colnamecheck(colnamevec, unigroupdf)

        if (isTRUE(all.equal(plottype, "percent"))) {
            gtypeaes <- ggplot2::aes(x = !!sym(kneename) / 2)
            if (is.na(xlimvec)) xlimvec <- c(0, 100)
            if (is.na(binwidthval)) binwidthval <- 5
            gtypelabs <- ggplot2::xlab("Distance TSS to knee (% of the gene)")
        } else {
            gtypeaes <- ggplot2::aes(x = (!!sym(kneename) * window_size) / 1000)
            if (is.na(xlimvec)) xlimvec <- c(0, 350)
            if (is.na(binwidthval)) binwidthval <- 10
            gtypelabs <- ggplot2::labs(x = "Distance TSS to knee (kb)",
                y = "Count",
                title = "TSS to knee position in kb for attenuated genes")
        }

        g <- ggplot2::ggplot(unigroupdf %>%
            dplyr::filter(!!sym(universename) == TRUE &
                !!sym(groupname) == "Attenuated"), gtypeaes) +
            ggplot2::geom_histogram(binwidth = binwidthval, fill = "grey",
                color = "black", boundary = 0) + ggplot2::xlim(xlimvec) +
            gtypelabs + ggplot2::theme_classic()

       if (plot) {
            warning("You chose to plot the auc, the figure is not saved.") # nolint
            print(g)
       } else {
            outfile <- paste0("histo_", plottype)
            ggplot2::ggsave(filename = paste0(outfile, ".", formatname),
                plot = g, device = formatname, path = outfold)
        }
}



##################
# MAIN
##################

## Reading objects and files
dfmeandiff <- readRDS(dfmeandiffpath)
unigroupdf <- readRDS(unigroupdfpath)
expdf <- read.csv(expdfpath, header = TRUE)


####
#### plotecdf
####

plotecdf(dfmeandiff, unigroupdf, expdf, "EGFR", colvec, outfold, plot = TRUE)
plotecdf(dfmeandiff, unigroupdf, expdf, "MARCHF6", colvec, outfold, plot = TRUE)


####
#### plotauc
####

genevec <- c("EGFR", "DAP", "FLI1", "MARCHF6", "LINC01619")
plotauc(unigroupdf, genevec, plot = TRUE)
plotauc(unigroupdf, legendpos = "none", subtitle = "Genes selected for Unibind",
    maintitle = "AUC Control vs HS", plot = TRUE, plottype = "groups")


####
#### metagene
####

plotmetagenes(unigroupdf, "attenuation", plot = TRUE)
plotmetagenes(unigroupdf, "outgroup", plot = TRUE)
plotmetagenes(unigroupdf, "universe", plot = TRUE)
plotmetagenes(unigroupdf, "all", plot = TRUE)


####
#### histogram
####

plothistoknee(unigroupdf, plottype = "percent", plot = TRUE)
plothistoknee(unigroupdf, plottype = "kb", plot = TRUE)
