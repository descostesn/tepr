library(dplyr)
library(purrr)
library(tidyr)



##################
# PARAMETERS
##################

working_directory <- "/g/romebioinfo/Projects/tepr/downloads"
extension <- "*.bg"
name_table <- "Cugusi_protein-lncRNA_stranded_analysis_MAPQ255_20230705.chr22.tsv" # nolint
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



main_table_read <- function(name_table, extension, working_directory,
    expression_threshold) {

    df_final <- data.frame()
    df_final_gene <- data.frame()
    concat_df <- data.frame()
    gene_name_list <- character()
    expressed_gene_name_list <- character()

    #Load data without header
    res <- getting_var_names(extension, working_directory)
    col_names <- res$col_names
    var_names <- res$var_names

    main_table <- data.frame()
    main_table <- read.delim(paste0(working_directory,"/", name_table),
                      header = FALSE, sep = "\t", col.names = col_names)

    for (gene in unique(main_table$gene)){
    gene_name_list <- c(gene_name_list, gene)
    }
    sorted_list <- sort(gene_name_list, decreasing = F)

    # Get the column names with the suffix "score"
    score_columns <- grep("score$", col_names, value = TRUE)

    # Initialize an empty list to store the mean column names
    mean_column_names <- list()
    expressed_gene_name <- data.frame()
    expressed_plus <- data.frame()

    expressed_transcript_name <- main_table %>%
    group_by(transcript) %>%
    dplyr::summarize(gene=gene[1],strand=strand[1],
                    across(all_of(score_columns), ~ mean(., na.rm = TRUE), .names = "{.col}_mean")) # nolint

    expressed_plus <- expressed_transcript_name %>%
    filter(strand=="+") %>% 
    select(gene, transcript, strand, contains("plus"))  %>%
    filter(across(all_of(contains("score")), ~ !is.na(.))) %>%
    filter(across(all_of(contains("score")), ~ . > expression_threshold))

    expressed_minus <- expressed_transcript_name %>%
    filter(strand == "-") %>% 
    select(gene, transcript, strand, contains("minus")) %>%
    filter(across(all_of(contains("score")), ~ !is.na(.))) %>%
    filter(across(all_of(contains("score")), ~ . > expression_threshold))

    expressed_transcript_name_list <- bind_rows(expressed_plus, expressed_minus) %>% arrange(transcript) %>% pull(transcript) # nolint

    return(list(main_table=main_table,expressed_transcript_name_list=expressed_transcript_name_list)) # nolint
}


genesECDF <- function(main_table, rounding, expressed_transcript_name_list,
    extension, working_directory) {

    gc()

    res <- getting_var_names(extension, working_directory)
    col_names <- res$col_names
    var_names <- res$var_names  

    total_iterations <- length(expressed_transcript_name_list)
    #setting the progress bar
    pb <- txtProgressBar(min = 0, max = total_iterations, style =5)

    j = 0
    concat_df <- data.frame()
    ## Looping through all the transcripts
    for (variable in expressed_transcript_name_list) {

        gene_table <- data.frame()
        bigDF <- data.frame()

        transcript_table <- data.frame()
        transcript <- filter(main_table, main_table$transcript == variable)
        transcript[transcript == "NAN"] <- NA
        bigDF <- transcript
        my_length <- length(bigDF[,'window'])

        var_names_score <- paste0(var_names,"_score")

        if (transcript$strand[1] == "-") {
            bigDF <- bigDF %>%
            select(!matches("plus"))
            bigDF$coord <- seq(from = my_length, to = 1, by = -1)
            bigDF <- arrange(bigDF, coord)
            conditions <- var_names_score[grepl("minus", var_names_score)]
        } else {
            bigDF <- bigDF %>%
            select(!matches("minus"))
            bigDF$coord <- seq(from=1, to=my_length,by=1)
            conditions <- var_names_score[grepl("plus",var_names_score)]
        }

        bigDF <- bigDF %>% fill(contains("score"), .direction = "downup")
        df_long <- bigDF %>% 
            gather(key = "variable", value = "value", conditions)
        df_long[,'value'] <- as.numeric(df_long[,'value'])
        df_long[,'value_round']<- round(df_long$value*rounding)

        j = j + 1
        #  Update the progress bar
        setTxtProgressBar(pb, j)

        list_df <- list()
        i <- 1
        for (my_var in unique(df_long$variable)) {
            df_subset <- subset(df_long, subset = variable == my_var) 
            df_expanded <- df_subset[rep(seq_len(nrow(df_subset)), df_subset$value_round), ] # nolint
            ecdf_df <- ecdf(df_expanded[,"coord"])
            df_subset$Fx <- ecdf_df(df_subset$coord) 
            list_df[[i]] <- df_subset
            i <- i + 1
        }

        df_final <- bind_rows(list_df)
        transcript_table <- df_final  %>% pivot_wider(., names_from = "variable", values_from = c("value", "value_round", "Fx")) %>% select(., -contains("value_round")) # nolint

        # getting rid of plus and minus
        if (transcript_table$strand[1]=="-") {
            # Drop columns containing "minus"
            columns_to_drop <- grep("plus", names(col_names), value = TRUE)
            dataset_without_dropped <- transcript_table %>%
            select(-all_of(columns_to_drop))

        # Modify column names by removing "_plus"
        modified_dataset <- dataset_without_dropped %>%
        rename_with(~gsub(".minus", "", .), contains(".minus"))
        } else {
            # Drop columns containing "minus"
            columns_to_drop <- grep("minus", names(col_names), value = TRUE)
            dataset_without_dropped <- transcript_table %>%
            select(-all_of(columns_to_drop))

            # Modify column names by removing "_plus"
            modified_dataset <- dataset_without_dropped %>%
            rename_with(~gsub(".plus", "", .), contains(".plus"))
        }
        concat_df <- bind_rows(concat_df, modified_dataset)
    }

    # # Close the progress bar
    close(pb)
    # list_gene_table <- concat_df %>% select(gene) %>% distinct()
    gc()

    return(concat_df = concat_df)
}



##################
# MAIN
##################

results_main_table <- main_table_read(name_table, extension, working_directory, 0.1) # nolint
resultsECDF <- genesECDF(main_table = results_main_table[[1]], rounding,
    expressed_transcript_name_list = results_main_table[[2]], extension,
    working_directory)
