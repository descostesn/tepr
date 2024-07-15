####################
# This script aims at ensuring the conformity of the R code with the bash one.
#
# Descostes - July 2024 - R-4.4.1
####################




##################
# PARAMETERS
##################

## Files obtained with bash
protcodbedshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/MANE_Select.protein_coding.bed" # nolint
lncrnabedshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/Ensembl_canonical_TSL123.lncRNA.bed" # nolint
protcodbednoblackwindshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/makewindow/v43.MANE_protein.window200.bed" # nolint
lncrnanednoblackwindshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/makewindow/v43.Ensembl_canonical_TSL123.lncRNA.bed" # nolint
blacklistshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/hg38-blacklist.v2.bed" # nolint
protcodnoblackfromshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/tmp2.bed" # nolint
lncrnanoblackfromshpath <- "/g/romebioinfo/Projects/tepr/downloads/annotations/tmp4.bed" # nolint
