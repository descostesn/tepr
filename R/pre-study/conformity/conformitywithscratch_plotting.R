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
colvec <- c("#90AFBB", "#FF9A04", "#10AFBB", "#FC4E07")
outfold <- "/g/romebioinfo/tmp/comparewithscratch-plotting"


##################
#FUNCTIONS
##################

.subtextvline <- function(condvec, geneinfo, digits) {

    reslist <- lapply(condvec, function(cond, geneinfo, digits) {

            ksname <- paste0("adjFDR_pvalaucks_", cond)
            ksval <- round(geneinfo[, ksname], digits)
            aucval <- round(geneinfo[, paste0("auc_", cond)], digits)

            if (ksval < 0.01) {
                kneeaucval <- geneinfo[, paste0("knee_AUC_", cond)]
                kneeval <- round(kneeaucval * geneinfo$windsize / 1000, digits)
                vlinecond <- data.frame(condition = cond, kneeval = kneeval)
            } else {
                kneeval <- "NA"
                vlinecond <- data.frame(condition = character(),
                    kneeval = numeric())
            }
            condtext <- paste0(cond, ": AUC = ", aucval, ", KS = ", ksval,
                ", Knee (kb) = ", kneeval)
            return(list(condtext, vlinecond))
        }, geneinfo, digits)
    subtext <- paste(sapply(reslist, "[", 1), collapse = "\n")
    vlinedf <- do.call("rbind", sapply(reslist, "[", 2))
    return(list(subtext, vlinedf))
}

.valcolbuild <- function(condvec, repvec) {
    return(as.vector(sapply(condvec,
        function(cond) { sapply(repvec, function(rep) {
            paste0(cond, rep, "score") })})))
}

.callggplot <- function(dflongecdf, colvec, windsizefact, vlinedf, subtext,
    outfold) {

    colvec <- as.vector(factor(dflongecdf$conditions, labels = colvec))
    ylimval <- 2 * max(dflongecdf$value)
    linexvals <- dflongecdf$coord * windsizefact
    lineyvals <- dflongecdf$coord / max(dflongecdf$coord)
    areayvals <- dflongecdf$value / ylimval

    g <- ggplot2::ggplot(dflongecdf, aes(x = coord, y = Fx,
        color = conditions)) +
        ggplot2::geom_line(aes(x = linexvals, y = lineyvals),
            linetype = "dashed", color = "red")

    g1 <- g + ggplot2::geom_area(aes(x = linexvals, y = areayvals,
        fill = conditions),
        alpha = 0.1,linewidth=0.2, position = 'identity') +
        scale_fill_manual(values = colvec) +
        geom_line(linewidth = 1, aes(x = linexvals)) +
        scale_color_manual(values = colvec)

    g2 <- g1 + ggplot2::scale_y_continuous(sec.axis = sec_axis(~ . * ylimval,
        name = "Transcription level")) +
        labs(x = "Distance from TSS (kb)",
             y = "Cumulative transcription density", title = genename,
             subtitle = subtext) + theme_classic()

    if (!isTRUE(all.equal(nrow(vlinedf), 0)))
        g2 <- g2 + ggplot2::geom_vline(data = vlinedf,
            aes(xintercept = kneeval),
            linetype = "dashed", color = "darkgrey")

    ggplot2::ggsave(filename = paste0(genename, ".pdf"),
        plot = g2, device = "pdf", path = outfold)
}

plotecdf <- function(dfmeandiff, unigroupdf, genename, colvec, outfold, # nolint
    digits = 2, middlewind = 100, verbose = TRUE) {

    ## Retrieving rows concerning the gene of interest
    if (verbose) message("\t Retrieving rows concerning the gene of interest")
    idxgene <- which(dfmeandiff$gene == genename)
    if (isTRUE(all.equal(length(idxgene), 0)))
        stop("The gene ", genename, " was not found")
    df <- dfmeandiff[idxgene, ]
    idxinfo <- which(unigroupdf$gene == genename)
    if (isTRUE(all.equal(length(idxinfo), 0)))
        stop("The gene ", genename, " was not found in unigroupdf")
    geneinfo <- unigroupdf[idxinfo, ]

    if (verbose) message("\t Gathering statistics about each condition")
    ## Computing the window size factor
    idxmiddle <- which(df$coord == middlewind)
    if (isTRUE(all.equal(length(idxmiddle), 0)))
        stop("The window ", middlewind, " was not found. Did you define a ",
        "number of windows equal to ", middlewind*2, "? If not, adjust the ",
        "parameter 'middlewind' to the half of the number of windows.")
    dfmid <- df[idxmiddle, ]
    windsizefact <- (dfmid$end - dfmid$start) / 1000

    ## Retrieving auc, ks, and knee
    condvec <- unique(expdf$condition)
    restmp <- .subtextvline(condvec, geneinfo, digits)
    subtext <- restmp[[1]]
    vlinedf <- restmp[[2]]

    ## Building data.frame for plotting with fx and value
    if (verbose) message("\t Building df for plotting with fx and value")
    repvec <- unique(expdf$replicate)
    colnamedfvec <- colnames(df)
    fxcolvec <- colnamedfvec[grep("^Fx_", colnamedfvec)]
    valcolvec <- .valcolbuild(condvec, repvec)
    ## Apply pivot
    dflongfx <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(fxcolvec), names_to = "conditions",
        values_to = "Fx") %>%
        dplyr::mutate(conditions = gsub("Fx_|_score", "", conditions)) # nolint
    dflongval <- df %>% tidyr::pivot_longer(
        cols = tidyselect::all_of(valcolvec), names_to = "conditions",
        values_to = "value") %>%
        dplyr::mutate(conditions = gsub("value_|_score", "", conditions)) # nolint
    ## merging
    commoncols <- intersect(names(dflongfx), names(dflongval))
    dflongecdf <- merge(dflongfx, dflongval, by = commoncols)

    ## Plotting
    if (verbose) message("\t Generating result to ", outfold)
    .callggplot(dflongecdf, colvec, windsizefact, vlinedf, subtext, outfold)
}



##################
# MAIN
##################

## Reading objects and files
dfmeandiff <- readRDS(dfmeandiffpath)
unigroupdf <- readRDS(unigroupdfpath)

####
#### plotecdf
####

plotecdf(dfmeandiff, unigroupdf, "EGFR", colvec, outfold)