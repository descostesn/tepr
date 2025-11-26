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
resecdflist <- genesECDF(resallexprs, expdf)

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
resecdflist <- genesECDF(resallexprs, expdf)

nbwindows <- resecdflist[[2]]
resecdf <- resecdflist[[1]]

saveRDS(resecdflist, file = "/Users/S238924/Documents/Other/PhD/PhD_papers/resecdflist_depleted_t0.1.dTAG") #threshold 0.1 for depleted ctrl and depleted hs


########################
## tepr and teprmulti

## Read input tables for WT 
expdf <- read.csv("expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE)
expdf$path=gsub("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/","",expdf$path)
transdf <- read.delim("dtag.tsv", header = FALSE)

## Calculate Average Expression and Filter Transcript Data
expthres <- 0.1

## Keep first
expdfone <- expdf %>% filter(condition=="wtctrl" | condition=="wtHS")
transdfone <- transdf %>% select(-(V10:V33))

## Keep second
expdftwo <- expdf %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdftwo <- transdf %>% select(-(V34:V57))
