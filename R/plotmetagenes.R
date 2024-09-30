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

plotmetagenes <- function(unigroupdf, dfmeandiff, plottype = "attenuation",
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
        ggplot2::geom_line(data = df, ggplot2::aes(x = coord / 2, # nolint
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
