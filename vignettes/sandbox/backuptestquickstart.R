library(tepr)

## For limiting
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpathall <- "/g/romebioinfo/tmp/preprocessing/objects-tsv-7cpus/cugusi.tsv"

## For testing
exppath <-  system.file("extdata", "exptab.csv", package="tepr")
transpath100 <- system.file("extdata", "cugusi_100.tsv", package="tepr")

expdf <- read.csv(exppath)
expthres <- 0.1

## For limiting
transdfall <- read.delim(transpathall, header = FALSE)

## For testing
transdf100 <- read.delim(transpath100, header = FALSE)
reslist100 <- tepr(expdf, transdf100, expthres, showtime = TRUE, verbose = TRUE)

## Limiting from all to two transcripts
genevec <- c("ENST00000275493.7", "ENST00000230895.11")
reslist <- lapply(genevec, function(currentgene, transdf){
    return(transdf[which(transdf$V5 == currentgene), ])
}, transdfall)
transdflim <- do.call("rbind", reslist)
reslistlim <- tepr(expdf, transdflim, expthres, showtime = TRUE, verbose = TRUE)

## Limiting from all to 6 transcripts
genevec <- c("ENST00000275493.7", "ENST00000230895.11", "ENST00000527786.7",
"ENST00000274140.10", "ENST00000549802.5", "ENST00000615891.2")
reslist <- lapply(genevec, function(currentgene, transdf){
    return(transdf[which(transdf$V5 == currentgene), ])
}, transdfall)
transdflim6 <- do.call("rbind", reslist)
reslistlim6 <- tepr(expdf, transdflim6, expthres, showtime = TRUE, verbose = TRUE)

## To enter tepr
alldf = transdflim
alldf = transdf100
alldf = transdflim6

expthres = 0.1; nbcpu = 1; rounding = 10; controlcondname = "ctrl"; stresscondname = "HS"; replaceval = NA
pval = 0.1; significant = FALSE; windsizethres = 50; countnathres = 20; meanctrlthres = 0.5; meanstressthres = 0.5
pvaltheorythres = 0.1; aucctrlthreshigher = -10; aucctrlthreslower = 15
aucstressthres = 15; attenuatedpvalksthres = 2; outgrouppvalksthres = 0.2; showtime = FALSE; verbose = TRUE

## To enter dfstrandlist <- mapply
strandname = unique(dfbytranscript$strand)[1]
directname = unique(expdf$strand)[1]
dfbytrans = dfbytranscript



genevec <- c("ENST00000275493.7", "ENST00000230895.11", "ENST00000527786.7",
"ENST00000274140.10", "ENST00000549802.5", "ENST00000615891.2")
genevec <- c("ENST00000275493.7", "ENST00000230895.11", "ENST00000527786.7",
"ENST00000274140.10", "ENST00000549802.5")
genevec <- c("ENST00000275493.7", "ENST00000230895.11", "ENST00000527786.7",
"ENST00000274140.10")
genevec <- c("ENST00000275493.7", "ENST00000230895.11", "ENST00000527786.7")


## To enter preprocessing
blacklistshpath = blacklistpath; nbcputrans = 1; finaltabpath = getwd()
finaltabname = "anno.tsv"; tmpfold = file.path(getwd(), "tmptepr")
saveobjectpath = NA; savefinaltable = TRUE; reload = FALSE
showtime = FALSE; showmemory = FALSE; deletetmp = TRUE; chromtab = NA
verbose = TRUE