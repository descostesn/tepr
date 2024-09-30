

# !!!!!!!!!!!!!!!
# write_file_protein_coding <- paste0(main_directory,"/protCoding_dTAG_Cugusi_stranded_20230810.tsv")
# write_file_lncRNA <- paste0(main_directory,"/lncRNA_dTAG_Cugusi_stranded_20230810.tsv")
# Big_tsv <-paste0(main_directory,"/dTAG_Cugusi_stranded_20230810.tsv")  
# main_directory <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/Xiong2021/analysis_KD/bedgraph"
# write_file_protein_coding <- paste0(main_directory,"/protCoding_XiongsiSAFB_stranded_20240214.tsv")
# write_file_lncRNA <- paste0(main_directory,"/lncRNA_XiongsiSAFB_stranded_20240214.tsv")
# Big_tsv <-paste0(main_directory,"/XiongsiSAFB_stranded_20240214.tsv")  
# main_directory <- "/Volumes/cristo-nas-shared/Personal_documents/Victor/DATA/SousaLuis2021/FASTQ_FILES/folder_clean/BAM/deduplicated/MAPQ255/stranded/bedgraph255"
# write_file_protein_coding <- paste0(main_directory,"/protCoding_SousaLuis_stranded_20240322.tsv")
# write_file_lncRNA <- paste0(main_directory,"/lncRNA_SousaLuis_stranded_20240322.tsv")
# Big_tsv <-paste0(main_directory,"/SousaLuis_stranded_20240322.tsv")  
# !!!!!!!!!!!!!!!!!!!!!!

library(dplyr)
library(purrr)

##################
# PARAMETERS
##################

main_directory <- "/g/romebioinfo/Projects/tepr/"
window <- 200
robjoutputfold <- "/g/romebioinfo/Projects/tepr/robjsave"


##################
#FUNCTIONS
##################

returnprotcodscoresfiles <- function(main_directory, subdirectory) {
    working_directory <- paste0(main_directory,"downloads/bedgraphs")
    bedgraph_files <- list.files(working_directory, pattern = "*.bg",
        full.names = TRUE)
    files <- bedgraph_files %>% map(~{
        filename <- tools::file_path_sans_ext(basename(.))
        file.path(working_directory, subdirectory,
        paste0(filename, ".window", window, ".MANE.wmean.name.score"))
    })
    return(files)
}

retrievescores <- function(bedgraph_files) {
    list_of_dfs <- lapply(bedgraph_files, read.delim, header = FALSE,
        sep = "\t", na.strings = "NAN", dec = ".",
        col.names = c("biotype", "chr", "coor1", "coor2", "transcript",
        "gene", "strand", "window", "id", "dataset", "score"),
        stringsAsFactors = FALSE)
    return(list_of_dfs)
}

##################
# MAIN
##################

#### protein coding

# reading all the files
bedgraph_files <- returnprotcodscoresfiles(main_directory,
    "protein_coding_score")
list_of_dfs <- retrievescores(bedgraph_files)

# joining all the files
# the last filter remove the PAR genes (pseudoautosomal genes both in X and Y)
joined_df <- purrr::reduce(list_of_dfs, dplyr::left_join,
    by = c("biotype", "chr", "coor1", "coor2", "transcript", "gene", "strand",
    "window", "id")) %>% dplyr::filter(strand != "Y")
saveRDS(joined_df, file = file.path(robjoutputfold, "joined_df.rds"))
write.table(joined_df,
 file = paste0(main_directory,
    "downloads/annotations/protCoding_dTAG_Cugusi_stranded_20230810.tsv"),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

rm(list_of_dfs)
rm(joined_df)

#### lncRNA

# reading all the files
bedgraph_files <- returnprotcodscoresfiles(main_directory,
    "lncRNA_score")
list_of_dfs <- retrievescores(bedgraph_files)

# joining all the files
joined_df_lncrna <- purrr::reduce(list_of_dfs, dplyr::left_join,
by = c("biotype","chr", "coor1", "coor2", "transcript", "gene", "strand",
"window", "id")) %>% dplyr::filter(strand != "Y")
saveRDS(joined_df_lncrna, file = file.path(robjoutputfold, "joined_df_lncrna.rds")) # nolint
write.table(joined_df_lncrna,
    file = paste0(main_directory,
        "downloads/annotations/lncRNA_dTAG_Cugusi_stranded_20230810.tsv"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)
rm(list_of_dfs)
rm(joined_df_lncrna)


### Merging the two files

working_directory <- paste0(main_directory, "downloads/annotations")
tsv_files <- list.files(working_directory, pattern = "*.tsv", full.names = TRUE)
files <- tsv_files %>% map(~{
    filename <- tools::file_path_sans_ext(basename(.))
    file.path(working_directory, paste0(filename, ".tsv"))
  })

# reading all the files
list_of_dfs <- lapply(files, read.delim, header = FALSE, sep = "\t",
    na.strings = "NAN", dec = ".", stringsAsFactors = FALSE)

# joining all the files
bound_df <- purrr::reduce(list_of_dfs, dplyr::bind_rows)
saveRDS(bound_df, file = file.path(robjoutputfold, "bound_df.rds"))

# write joined table
big_tsv <- paste0(working_directory, "/dTAG_Cugusi_stranded_20230810.tsv")
write.table(bound_df, file = big_tsv, sep = "\t", row.names = FALSE,
    col.names = FALSE, quote = FALSE)

rm(list_of_dfs)
rm(bound_df)
