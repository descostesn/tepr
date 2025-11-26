########################
## FOR WT

## Read input tables for WT 
expdf <- read.csv("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE) %>% filter(condition=="wtctrl" | condition=="wtHS")
transdf <- read.delim("/Users/S238924/Documents/Other/PhD/PhD_papers/dtag.tsv", header = FALSE) %>% select(-(V10:V33))

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

expdf <- read.csv("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE) %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdf <- read.delim("/Users/S238924/Documents/Other/PhD/PhD_papers/dtag.tsv", header = FALSE) %>% select(-(V34:V57))

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

