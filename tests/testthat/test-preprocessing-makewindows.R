## Data
exptabpath <- system.file("extdata", "exptab-preprocessing.csv", package="tepr")
gencodepath <- system.file("extdata", "gencode-chr13.gtf", package = "tepr")
windsize <- 200

## Necessary result to call makewindows
allannobed <- retrieveanno(exptabpath, gencodepath)

!!!!!!!!!!

allwindowsbed <- makewindows(allannobed, windsize, nbcputrans,
                verbose, saveobjectpath, showtime)

