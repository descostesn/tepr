########################
## FOR WT

## Read input tables for WT 
expdf <- read.csv("expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE) %>% filter(condition=="wtctrl" | condition=="wtHS")
expdf$path=gsub("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/","",expdf$path)

transdf <- read.delim("dtag.tsv", header = FALSE) %>% select(-(V10:V33))

## Calculate Average Expression and Filter Transcript Data
expthres <- 0.1
resallexprs <- averageandfilterexprs(expdf, transdf, expthres)

## Count NA values per transcript and condition
rescountna <- countna(resallexprs, expdf)
## Compute ECDF for Genes Based on Expression Data
resecdflist <- genesECDF(resallexprs)

nbwindows <- resecdflist[[2]]
resecdf <- resecdflist[[1]]

saveRDS(resecdflist, file = "/Users/S238924/Documents/Other/PhD/PhD_papers/resecdflist_wt_t0.1.dTAG") #threshold 0.1 for wt ctrl and wt hs

####################################
## Read input tables for DEPLETED ##
###############################

expdf <- read.csv("expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE) %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdf <- read.delim("dtag.tsv", header = FALSE) %>% select(-(V34:V57))
expdf$path=gsub("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/","",expdf$path)

## Calculate Average Expression and Filter Transcript Data
expthres <- 0.1
resallexprs <- averageandfilterexprs(expdf, transdf, expthres)

## Count NA values per transcript and condition
rescountna <- countna(resallexprs, expdf)
## Compute ECDF for Genes Based on Expression Data
resecdflist <- genesECDF(resallexprs)

nbwindows <- resecdflist[[2]]
resecdf <- resecdflist[[1]]

saveRDS(resecdflist, file = "/Users/S238924/Documents/Other/PhD/PhD_papers/resecdflist_depleted_t0.1.dTAG") #threshold 0.1 for depleted ctrl and depleted hs


########################
## Testing tepr and teprmulti
########################

## Read input tables for WT 
expdf <- read.csv("expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE)
expdf$path=gsub("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/","",expdf$path)
transdf <- read.delim("dtag.tsv", header = FALSE)

## Calculate Average Expression and Filter Transcript Data
expthres <- 0.1

###
# Tepr on both conditions separately
###

## Keep first
expdfone <- expdf %>% filter(condition=="wtctrl" | condition=="wtHS")
transdfone <- transdf %>% select(-(V10:V33))
resone <- tepr(expdfone, alldf=transdfone, expthres, controlcondname = "wtctrl", stresscondname = "wtHS")

## Keep second
expdftwo <- expdf %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdftwo <- transdf %>% select(-(V34:V57))
restwo <- tepr(expdftwo, alldf=transdftwo, expthres, controlcondname = "depletedctrl", stresscondname = "depletedHS")

saveRDS(expdftwo, file="expdftwo.RDS")
saveRDS(transdftwo, file="transdftwo.RDS")

###
# teprmulti on all conditions
###

## teprmulti (error occurs here)
resmulti <- teprmulti(expdf, alldf=transdf, expthres)

## Debug

alldf=transdf; nbcpu = 1; rounding = 10; dontcompare = NULL; replaceval = NA; pval = 0.1; significant = FALSE; windsizethres = 50; countnathres = 20; pvaltheorythres = 0.1; meancond1thres = 0.5; meancond2thres = 0.5; auccond1threshigher = -10; auccond1threslower = 15; auccond2thres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2; saveobjectpath = NA; reload = FALSE; showtime = FALSE; showmemory = FALSE; verbose = TRUE

dontcompare <- c("depletedctrl_vs_wtctrl", "depletedctrl_vs_wtHS", "depletedHS_vs_wtctrl", "depletedHS_vs_wtHS")

currentcol <- matcond[ ,1]



tepr(expdf = expdf2cond, alldf = alldf2cond, expthres = expthres, nbcpu = nbcpu, rounding = rounding,
            controlcondname = cond1name, stresscondname = cond2name,
            replaceval = replaceval, pval = pval, significant = significant,
            windsizethres = windsizethres, countnathres = countnathres,
            meanctrlthres = meancond1thres, meanstressthres = meancond2thres,
            pvaltheorythres = pvaltheorythres,
            aucctrlthreshigher = auccond1threshigher,
            aucctrlthreslower = auccond1threslower,
            aucstressthres = auccond2thres,
            attenuatedpvalksthres = attenuatedpvalksthres,
            outgrouppvalksthres = outgrouppvalksthres, showtime = showtime,
            verbose = verbose)

tepr(expdf = expdf2cond, alldf = alldf2cond, expthres = expthres, controlcondname = cond1name, stresscondname = cond2name)

expdftwo <- expdf %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdftwo <- transdf %>% select(-(V34:V57))

isTRUE(all.equal(expdf2cond, expdftwo))
isTRUE(all.equal(alldf2cond, transdftwo))

# Set up data for tepr function testing
colnames(alldf2cond) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33")
