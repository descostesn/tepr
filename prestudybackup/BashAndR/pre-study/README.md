This script contains the code of the documentation (see corresponding folder).

It enabled to re-perform the pre-processing steps and to save files that were used to ensure the conformity of the R code with the bash one. See the script [preprocessing.R](../../R/pre-study/preprocessing.R) that needs to be run first. This script produces different R objects that are saved to the path defined by the variable `robjoutputfold`.

These objects are then used in the script [conformity_with_bash.R](../../R/pre-study/conformity_with_bash.R).

The script [test-preprocessing.R](test-preprocessing.R) evaluates the R code contained in [TSV_creation_lite_20240617.html](documentation/TSV_creation_lite_20240617.html).
