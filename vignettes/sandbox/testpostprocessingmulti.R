library("tepr")

##################
# PARAMETERS
##################

exptabpath <- "/g/romebioinfo/Projects/tepr-data/downloads/annotations/exptab-bedgraph-DRB.csv" # nolint
finaltabpath <- "/g/romebioinfo/tmp/preprocessing-drbseq/drbttseq.tsv"
finaltabpathvic <- "/g/romebioinfo/Projects/tepr-data/downloads/bedgraphs-DRB/DRB_Cugusi_20240619.tsv"
expthres <- 0.1
nbcpu <- 5
rounding <- 10
dontcompare <- NULL
controlcondname <- "ctrl"
stresscondname <- "HS"
replaceval <- NA
pval <- 0.1
significant <- FALSE
windsizethres <- 50
countnathres <- 20
meancond1thres <- 0.5
meancond2thres <- 0.5
pvaltheorythres <- 0.1
auccond1threshigher <- -10
auccond1threslower <- 15
auccond2thres <- 15
attenuatedpvalksthres <- 2
outgrouppvalksthres <- 0.2
showtime <- TRUE
verbose <- TRUE
saveobjectpath <- "/g/romebioinfo/tmp/multitest"

## matcond <- matcond[, c(1, 6, 10, 24)]



##################
# MAIN
##################

expdf <- read.csv(exptabpath, header = TRUE)
alldf <- read.delim(finaltabpath, header = FALSE)

expdf <- expdf[which(expdf$condition == "ctrl10" | expdf$condition == "ctrl20" |
 expdf$condition == "HS20"), ]
resteprmulti <- teprmulti(expdf, alldf, expthres, nbcpu = 5, showtime = TRUE,
    saveobjectpath = "/g/romebioinfo/tmp/testmultifun/objbackup")
# resteprmulti <- readRDS("/g/romebioinfo/tmp/testmultifun/resteprmulti.rds")


# digits = 2; middlewind = 100; pval = 0.01
# colvec = c("#90AFBB", "#10AFBB", "#FF9A04", "#FC4E07")
# genaucvec = NA; aucaxisminx = -10; aucaxismaxx = 100; aucaxisminy = -10
# aucaxismaxy = 100; aucmaintitle = ""; aucsubtitle = ""
# auclegendpos = "bottom"; formatname = "pdf"; uniname = "Universe"
# groupname = "Group"; histkneexlim = NA; binwidthvalhistknee = NA

plotmulti(resteprmulti, expdf, ecdfgenevec = c("CDC27", "BCAR1", "TRAM2"),
    outfold = "/g/romebioinfo/tmp/testmultifun")

## test plotecdf
dfmeandiff = complist[[1]]; unigroupdf = complist[[2]]
genename = currentgene; outfold = outfoldcomp; plot = FALSE

## test plotauc
tab = complist[[2]]; genevec = genaucvec;auc_ctrlname = name1
auc_stressname = name2;pvalkstestcolname = pvalks; axismin_x = aucaxisminx;
axismax_x = aucaxismaxx; axismin_y = aucaxisminy;axismax_y = aucaxismaxy
maintitle = aucmaintitle;subtitle = aucsubtitle; legendpos = auclegendpos
outfold = outfoldcomp;outfile = aucfilename; plottype = "groups"
plot = FALSE; universename = uniname; groupname = groupname
            



## vic input for figures
# AUC_knee_DRB_vic <- read.delim( file = "/g/romebioinfo/Projects/tepr-data/downloads/inputfiles-DRBanalysis/AUC_knee_DRB.tsv", sep = "\t",  header = TRUE)

# ## subset of the result of tepr multi
# reslist <- readRDS(file.path(saveobjectpath, "reslist.rds"))
# resteprmulti <- reslist

#alldfvic <- read.delim(finaltabpathvic, header = FALSE)

# expdf_backup = expdf
# alldf_backup = alldf
# alldf_backupvic = alldfvic
# expdf = expdf2cond

# alldf = alldf2cond
# alldfvic = alldf2condvic

# controlcondname = cond1name
# stresscondname = cond2name; meanctrlthres = meancond1thres
# meanstressthres = meancond2thres; aucctrlthreshigher = auccond1threshigher
# aucctrlthreslower = auccond1threslower; aucstressthres = auccond2thres
