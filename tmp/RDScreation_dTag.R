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
expthres <- 0.1

# Victor
resallexprs <- averageandfilterexprs(expdf, transdf, expthres)
## Count NA values per transcript and condition
rescountna <- countna(resallexprs, expdf)
## Compute ECDF for Genes Based on Expression Data
resecdflist <- genesECDF(resallexprs)
nbwindows <- resecdflist[[2]]
resecdf <- resecdflist[[1]]
saveRDS(resecdflist, file = "/Users/S238924/Documents/Other/PhD/PhD_papers/resecdflist_depleted_t0.1.dTAG") #threshold 0.1 for depleted ctrl and depleted hs

################################################# debug
#tepr
alldf=transdf;rounding = 10; dontcompare = NULL; replaceval = NA; pval = 0.1; significant = FALSE; windsizethres = 50; countnathres = 20; pvaltheorythres = 0.1; meancondonethres = 0.5; meancondtwothres = 0.5; auccondonethreshigher = -10; auccondonethreslower = 15; auccondtwothres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2; saveobjectpath = NA; reload = FALSE; showtime = FALSE; showmemory = FALSE; verbose = TRUE
# save resallexprs[[1]]
saveRDS(resallexprs[[1]], file="tworesallexprs1.rds")
saveRDS(resallexprs[[2]], file="tworesallexprs2.rds")


########################
## Testing tepr and teprmulti
########################

## Read input tables for WT 
expdf <- read.csv("expdf_dTag_corrected.csv", header = T, stringsAsFactors = FALSE)
expdf$path=gsub("/Users/S238924/Documents/Other/PhD/PhD_papers/dTAG/","",expdf$path)
transdf <- read.delim("dtag.tsv", header = FALSE)
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
dontcompare <- c("depletedctrl_vs_wtctrl", "depletedctrl_vs_wtHS", "depletedHS_vs_wtctrl", "depletedHS_vs_wtHS")
resmulti <- teprmulti(expdf, alldf=transdf, expthres, dontcompare = dontcompare)

######################################################################## Debug
# teprmulti
alldf=transdf; nbcpu = 1; rounding = 10; dontcompare = NULL; replaceval = NA; pval = 0.1; significant = FALSE; windsizethres = 50; countnathres = 20; pvaltheorythres = 0.1; meancondonethres = 0.5; meancondtwothres = 0.5; auccondonethreshigher = -10; auccondonethreslower = 15; auccondtwothres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2; saveobjectpath = NA; reload = FALSE; showtime = FALSE; showmemory = FALSE; verbose = TRUE
dontcompare <- c("depletedctrl_vs_wtctrl", "depletedctrl_vs_wtHS", "depletedHS_vs_wtctrl", "depletedHS_vs_wtHS")
currentcol <- matcond[ ,1]
# tepr
expdf = expdftwocond; alldf = alldftwocond; expthres = expthres; nbcpu = nbcpu; rounding = rounding; controlcondname = condonename; stresscondname = condtwoname; replaceval = replaceval; pval = pval; significant = significant; windsizethres = windsizethres; countnathres = countnathres; meanctrlthres = meancondonethres; meanstressthres = meancondtwothres; pvaltheorythres = pvaltheorythres; aucctrlthreshigher = auccondonethreshigher; aucctrlthreslower = auccondonethreslower; aucstressthres = auccondtwothres; attenuatedpvalksthres = attenuatedpvalksthres; outgrouppvalksthres = outgrouppvalksthres; showtime = showtime; verbose = verbose
# compare with manual two conditions
tworesallexprs1 <- readRDS("tworesallexprs1.rds")
tworesallexprs2 <- readRDS("tworesallexprs2.rds")

#TRUE
identical(resallexprs[[1]], tworesallexprs1)

!!# FALSE
!! identical(resallexprs[[2]], tworesallexprs2)














tepr(expdf = expdftwocond, alldf = alldftwocond, expthres = expthres, nbcpu = nbcpu, rounding = rounding,
            controlcondname = cond1name, stresscondname = cond2name,
            replaceval = replaceval, pval = pval, significant = significant,
            windsizethres = windsizethres, countnathres = countnathres,
            meanctrlthres = meancondonethres, meanstressthres = meancondtwothres,
            pvaltheorythres = pvaltheorythres,
            aucctrlthreshigher = auccondonethreshigher,
            aucctrlthreslower = auccondonethreslower,
            aucstressthres = auccondtwothres,
            attenuatedpvalksthres = attenuatedpvalksthres,
            outgrouppvalksthres = outgrouppvalksthres, showtime = showtime,
            verbose = verbose)

tepr(expdf = expdftwocond, alldf = alldftwocond, expthres = expthres, controlcondname = condonename, stresscondname = condtwoname)

expdftwo <- expdf %>% filter(condition=="depletedctrl" | condition=="depletedHS")
transdftwo <- transdf %>% select(-(V34:V57))

isTRUE(all.equal(expdftwocond, expdftwo))
isTRUE(all.equal(alldftwocond, transdftwo))

identical(expdftwocond, expdftwo)
identical(ncol(transdftwo), ncol(alldftwocond))
backupcolnames <-colnames(alldftwocond)
colnames(alldftwocond) <- colnames(transdftwo)
identical(alldftwocond, transdftwo)
identical("depletedctrl", condonename)
identical("depletedHS", condtwoname)

# Set up data for tepr function testing
colnames(alldftwocond) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33")
