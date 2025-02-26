exptabpath = "newfilesreduced/exptab_limited.csv"
gencodepath = "newfilesreduced/gencode-chr13_limited.gtf"
windsize = 200
maptrackpath = "k50.umap.chr13.hg38.0.8.bed"
blacklistpath = "hg38-blacklist-chr13.v2.bed"
genomename = "hg38"
nbcputrans = 5; finaltabpath = "./"; finaltabname = "anno.tsv"
tmpfold = "./tmp"; saveobjectpath = NA; savefinaltable = TRUE;
reload = FALSE; showtime = TRUE; showmemory = TRUE; deletetmp = TRUE
chromtab = NA; verbose = TRUE

## makewindows
expbed = allannobed; nbwindows = windsize

