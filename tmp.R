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
