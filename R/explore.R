library(dplyr)
library(purrr)


##################
# PARAMETERS
##################

working_directory <- "/g/romebioinfo/Projects/tepr/downloads" 
extension <- "*.bg"
name_table <- "Cugusi_protein-lncRNA_stranded_analysis_MAPQ255_20230705.chr22.tsv"
rounding <- 10


##################
#FUNCTIONS
##################


getting_var_names <- function(extension, working_directory) {

  bedgraph_files <- list.files(working_directory, pattern = extension,
    full.names = TRUE)
  files <- bedgraph_files %>% map(~{filename <- tools::file_path_sans_ext(basename(.))}) # nolint
  string <- files
  var_names <- string

  for (i in seq_along(var_names)) {
    if (grepl("(reverse|forward)", var_names[i])) {
      var_names[i] <- gsub("reverse", "minus", var_names[i])
      var_names[i] <- gsub("forward", "plus", var_names[i])
    }
  }

  # Extract conditions
  Conditions <- unique(sub("(\\w+)_rep\\d+.*", "\\1", var_names)) ## verify it can work with several "_" # nolint

  # Extract replication numbers
  replicate_numbers <- unique(sub(".*_rep(\\d+).*", "\\1", var_names))

  # Setting fixed column names
  fixed_names <- c("biotype","chr", "coor1", "coor2","transcript", "gene", "strand","window","id") #biotype is added in the .tsv file # nolint

  num_fixed_cols <- length(fixed_names)
  # Number of variable columns based on input file
  num_var_cols <- ncol(table) - num_fixed_cols
  #Generate variable column names for category 2, alternating with category 1
  test <- rep(var_names, each=2)
  suffix <- rep(c("", "_score"), length(var_names))
  test <- paste0(test, suffix)
  col_names <- c(fixed_names, test)

  return(list(col_names=col_names,var_names=var_names, replicate_numbers=replicate_numbers, Conditions=Conditions)) # nolint
}


##################
# MAIN
##################

getting_var_names(extension, working_directory)
