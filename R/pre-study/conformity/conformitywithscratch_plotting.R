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
