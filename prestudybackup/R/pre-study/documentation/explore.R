library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

# !!!!!!!
# environment: /g/romebioinfo/tmp/explore
# !!!!!!!!!


##################
# PARAMETERS
##################

#working_directory <- "/g/romebioinfo/Projects/tepr/downloads"
working_directory <- "/mnt/c/Users/descoste/Documents/analysis/cristofari/explore"
extension <- "*.bg"
name_table <- "/g/romebioinfo/Projects/tepr/downloads/annotations/dTAG_Cugusi_stranded_20230810.tsv" # nolint
rounding <- 10


##################
#FUNCTIONS
##################


getting_var_names <- function(extension, workdir) {

  bedgraph_files <- list.files(workdir, pattern = extension,
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




main_table_read <- function(name_table, extension, workingdirectory,
    expression_threshold) {

    df_final <- data.frame()
    df_final_gene <- data.frame()
    concat_df <- data.frame()
    gene_name_list <- character()
    expressed_gene_name_list <- character()

    #Load data without header
    res <- getting_var_names(extension, workingdirectory)
    col_names <- res$col_names
    var_names <- res$var_names

    main_table <- data.frame()
    main_table <- read.delim(name_table, header = FALSE, sep = "\t",
      col.names = col_names)

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

genesECDF <- function(main_table, rounding, expressed_transcript_name_list, extension, working_directory){
gc()
  
res <- getting_var_names(extension, working_directory)
col_names <- res$col_names
var_names <- res$var_names  
  
total_iterations <- length(expressed_transcript_name_list)
#setting the progress bar
pb <- txtProgressBar(min = 0, max = total_iterations, style =5)

j=0
concat_df <- data.frame()
## Looping through all the transcripts
for (variable in expressed_transcript_name_list){
  gene_table <- data.frame()
  bigDF <- data.frame()  
  
  transcript_table <- data.frame()
  transcript <- filter(main_table, main_table$transcript==variable) 
  transcript[transcript == "NAN"] <- NA
  bigDF <- transcript
  my_length <- length(bigDF[,'window'])
  
  
  var_names_score <- paste0(var_names,"_score")
  
  if (transcript$strand[1]=="-") {
    bigDF <- bigDF %>%
      select(!matches("plus"))
    bigDF$coord <- seq(from=my_length, to=1,by=-1)
    bigDF <- arrange(bigDF, coord)
    conditions <- var_names_score[grepl("minus",var_names_score)]
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
  
  j=j+1
#  Update the progress bar
  setTxtProgressBar(pb, j)
  
  list_df <- list()
  i <- 1
  for (my_var in unique(df_long$variable)){
    
    df_subset <- subset(df_long, subset = variable == my_var) 
    df_expanded <- df_subset[rep(seq_len(nrow(df_subset)), df_subset$value_round), ]
    ecdf_df <- ecdf(df_expanded[,"coord"])
    df_subset$Fx <- ecdf_df(df_subset$coord) 
    list_df[[i]] <- df_subset
    i <- i + 1
  }
  
  df_final <- bind_rows(list_df)
  transcript_table <- df_final  %>% pivot_wider(., names_from = "variable", values_from = c("value", "value_round", "Fx")) %>% select(., -contains("value_round")) 
  
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

return(concat_df=concat_df)
}

# genesECDF <- function(main_table, rounding, expressed_transcript_name_list,
#     extension, workdir) {

#     gc()

#     res <- getting_var_names(extension, workdir)
#     col_names <- res$col_names
#     var_names <- res$var_names  

#     total_iterations <- length(expressed_transcript_name_list)
#     #setting the progress bar
#     pb <- txtProgressBar(min = 0, max = total_iterations, style =5)

#     j = 0
#     concat_df <- data.frame()

#     ## Looping through all the transcripts i.e performs the loop for each
#     ## transcript that has been kept by the function main_table_read because
#     ## it was expressed
#     for (variable in expressed_transcript_name_list) {

#         gene_table <- data.frame()
#         bigDF <- data.frame()

#         ## Isolating the rows corresponding to the transcript
#         ## Question: Could we do it at once for all expressed transcripts with
#         ## a match (transcriptmaintable in selectedtrans, remove na, use the
#         ## transcript column as factor to split the main table)
#         transcript_table <- data.frame()
#         transcript <- filter(main_table, main_table$transcript == variable)
#         ## This replacement is not necessary in my hands
#         transcript[transcript == "NAN"] <- NA
#         ## Changing the name of transcript can be useful to conserve it. So far,
#         ## it seems to be done to keep looking for the strand. It could be
#         ## stored in a variable. Add verification that the strand is unique
#         ## in transcript.
#         bigDF <- transcript
#         ## This is equivalent to nrow(transcript). Check if every transcript has
#         ## the same number of windows, if it is the case, it is not necessary
#         ## to take this into account.
#         my_length <- length(bigDF[,'window'])

# # ---------------------------------------------------------------------------------------------------
# ## SUMMARY
# ##
# ## This section filters the score columns according to the strand of the transcript considered

#         ## This only builds the names of the columns containing the scores
#         var_names_score <- paste0(var_names,"_score")

#         if (transcript$strand[1] == "-") { # see above, the strand can be stored in a variable instead of calculating it. # nolint
#             bigDF <- bigDF %>%
#             select(!matches("plus"))
#             bigDF$coord <- seq(from = my_length, to = 1, by = -1)
#             bigDF <- arrange(bigDF, coord)
#             conditions <- var_names_score[grepl("minus", var_names_score)]
#         } else {
#             bigDF <- bigDF %>%
#             select(!matches("minus"))
#             bigDF$coord <- seq(from=1, to=my_length,by=1)
#             conditions <- var_names_score[grepl("plus",var_names_score)]
#         }
# # ---------------------------------------------------------------------------------------------------

#         ## To my understanding this Fills missing values in selected columns
#         ## using first down and then up previous entry.
#         bigDF <- bigDF %>% fill(contains("score"), .direction = "downup")

#         ## The code below creates a data.frame of two columns with the name of
#         ## the experiment as key and the expression as value. This is the kind
#         ## of transformation we do before using ggplot2. However I do not see
#         ## why keeping the previous columns is necessary. They will be repeated
#         ## several times.
#         df_long <- bigDF %>% 
#             gather(key = "variable", value = "value", conditions)
#         ## Makes the values as integer by multiplying by 10. Why not using the
#         ## ceil or floor function. Is it important to increase the scores?
#         df_long[,'value'] <- as.numeric(df_long[,'value'])
#         df_long[,'value_round']<- round(df_long$value*rounding)

#         #  Update the progress bar
#         j = j + 1
#         setTxtProgressBar(pb, j)

#         list_df <- list()
#         i <- 1
#         ## This for loop is equivalent to computing on each column, bigDF might
#         ## not be useful.
#         for (my_var in unique(df_long$variable)) {
#             df_subset <- subset(df_long, subset = variable == my_var) # This is just selecting the lines that we had in the initial table # nolint
#             df_expanded <- df_subset[rep(seq_len(nrow(df_subset)), df_subset$value_round), ] # nolint
#             ecdf_df <- ecdf(df_expanded[,"coord"])
#             df_subset$Fx <- ecdf_df(df_subset$coord) 
#             list_df[[i]] <- df_subset
#             i <- i + 1
#         }

#         ## This two lines are equivalent to cbind a matrix of Fx to the transcript table used at the beginning
#         ## The columns are biotype, chr, coor1, coor2, transcript, gene, strand, window, id, ctrl_rep1.plus, ctrl_rep2.plus,
#         ## HS_rep1.plus, HS_rep2.plus, coord, value_ctrl_rep1.plus_score, value_ctrl_rep2.plus_score,
#         ## value_HS_rep1.plus_score, value_HS_rep2.plus_score, Fx_ctrl_rep1.plus_score, Fx_ctrl_rep2.plus_score,
#         ## Fx_HS_rep1.plus_score, Fx_HS_rep2.plus_score
#         df_final <- bind_rows(list_df)
#         transcript_table <- df_final  %>% pivot_wider(., names_from = "variable", values_from = c("value", "value_round", "Fx")) %>% select(., -contains("value_round")) # nolint

# #---------------
# ## ADDED BY ME
# #modified_dataset <- transcript_table
# #---------------


# # ------------------------------------------------
# ## THIS CANNOT BE TRIGGERED BECAUSE THE COLUMNS DO NOT EXIST (????)
#         # getting rid of plus and minus
#         if (transcript_table$strand[1]=="-") {
#             # Drop columns containing "minus"
#             columns_to_drop <- grep("plus", col_names, value = TRUE)
#             dataset_without_dropped <- transcript_table %>%
#             select(-all_of(columns_to_drop))

#         # # Modify column names by removing "_plus"
#         modified_dataset <- dataset_without_dropped %>%
#         rename_with(~gsub(".minus", "", .), contains(".minus"))
#         } else {
#             # Drop columns containing "minus"
#             columns_to_drop <- grep("minus", col_names, value = TRUE)
#             dataset_without_dropped <- transcript_table %>%
#             select(-all_of(columns_to_drop))

#         #     # Modify column names by removing "_plus"
#             modified_dataset <- dataset_without_dropped %>%
#             rename_with(~gsub(".plus", "", .), contains(".plus"))
#          }
# # ----------------------------------------------------
#         concat_df <- bind_rows(concat_df, modified_dataset)
#     }
# #saveRDS(concat_df, file = file.path("/g/romebioinfo/Projects/tepr/robjsave/concatdf_fromexplore.rds"))
#     # # Close the progress bar
#     close(pb)
#     # list_gene_table <- concat_df %>% select(gene) %>% distinct()
#     gc()

#     return(concat_df = concat_df)
# }


# calculates_meanFx <- function(concat_df,window_number){

# res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
# Conditions <- res$Conditions
# replicate_numbers <- res$replicate_numbers


# column_vector_value <- character()
# column_vector_Fx <- character()


# for (cond in Conditions) {
#   mean_value_condi_name <- paste0("mean_value_", cond)
#   mean_Fx_condi_name <- paste0("mean_Fx_", cond)
#   diff_Fx_condi_name <- paste0("diff_Fx_", cond)
  
#   for (rep_num in replicate_numbers) {
#     new_column_value <- paste0("value_", cond, "_rep", rep_num, "_score") # Generate a new item
#     new_column_Fx <- paste0("Fx_", cond, "_rep", rep_num, "_score") # Generate a new item
#     column_vector_value <- c(column_vector_value, new_column_value)
#     column_vector_Fx <- c(column_vector_Fx, new_column_Fx)
#   }
#   # Calculate row means for the specified columns
#   # Check if there is more than one replicate
#   if (length(replicate_numbers) > 1) {
#     concat_df[[mean_value_condi_name]] <- rowMeans(concat_df[, column_vector_value], na.rm = F)
#     concat_df[[mean_Fx_condi_name]] <- rowMeans(concat_df[, column_vector_Fx], na.rm = FALSE)
#   } else {
#     # Handle case when column_vector_value is empty
#     new_column_value <- paste0("value_", cond, "_rep", "1", "_score") # Generate a new item
#     new_column_Fx <- paste0("Fx_", cond, "_rep", "1", "_score") # Generate a new item
    
#     concat_df[[mean_value_condi_name]] <- concat_df[[new_column_value]]
#     concat_df[[mean_Fx_condi_name]] <- concat_df[[new_column_Fx]]
#   }
  
#   concat_df[[diff_Fx_condi_name]] <- concat_df[[mean_Fx_condi_name]] - concat_df$coord/window_number ## Difference with the y=x ECDF, used to calculate AUC
  
#   column_vector_value <- character() ## obligatory to reset the columns values to empty
#   column_vector_Fx <- character()
# }

# return(concat_dfFx=concat_df)
# }

# concat_df <- resultsECDF
# window_number <- 200
calculates_meanFx <- function(concat_df,window_number){

res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
Conditions <- res$Conditions
replicate_numbers <- res$replicate_numbers


column_vector_value <- character()
column_vector_Fx <- character()


for (cond in Conditions) {
  mean_value_condi_name <- paste0("mean_value_", cond)
  mean_Fx_condi_name <- paste0("mean_Fx_", cond)
  diff_Fx_condi_name <- paste0("diff_Fx_", cond)
  
  for (rep_num in replicate_numbers) {
    new_column_value <- paste0("value_", cond, "_rep", rep_num, "_score") # Generate a new item
    new_column_Fx <- paste0("Fx_", cond, "_rep", rep_num, "_score") # Generate a new item
    column_vector_value <- c(column_vector_value, new_column_value)
    column_vector_Fx <- c(column_vector_Fx, new_column_Fx)
  }
  # Calculate row means for the specified columns
  # Check if there is more than one replicate
  if (length(replicate_numbers) > 1) {
    concat_df[[mean_value_condi_name]] <- rowMeans(concat_df[, column_vector_value], na.rm = F)
    concat_df[[mean_Fx_condi_name]] <- rowMeans(concat_df[, column_vector_Fx], na.rm = FALSE)
  } else {
    # Handle case when column_vector_value is empty
    new_column_value <- paste0("value_", cond, "_rep", "1", "_score") # Generate a new item
    new_column_Fx <- paste0("Fx_", cond, "_rep", "1", "_score") # Generate a new item
    
    concat_df[[mean_value_condi_name]] <- concat_df[[new_column_value]]
    concat_df[[mean_Fx_condi_name]] <- concat_df[[new_column_Fx]]
  }
  
  concat_df[[diff_Fx_condi_name]] <- concat_df[[mean_Fx_condi_name]] - concat_df$coord/window_number ## Difference with the y=x ECDF, used to calculate AUC
  
  column_vector_value <- character() ## obligatory to reset the columns values to empty
  column_vector_Fx <- character()
}

return(concat_dfFx=concat_df)
}



condition_comparison <- function(extension,working_directory) {

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

  for (i in 1:length(Conditions)) {
    for (j in (1+i):length(Conditions)) {
      if (j>length(Conditions)){ break}
      cond1 <- Conditions[i]
      cond2 <- Conditions[j]
      newcol <- paste0(cond1," vs ", cond2)
      newcol_reverse <- paste0(cond2," vs ", cond1)
      print(newcol)
 #   print(newcol_reverse)
    }
  }
}

condition_compared <- function(extension, working_directory,
  dontcompare = NULL) {

    res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
    Conditions <- res$Conditions

    for (i in 1:length(Conditions)) {
      cond1 <- Conditions[i]
      #  print(cond1)
      for (j in (1+i):length(Conditions)) {
        if (j>length(Conditions)){ break}
        cond2 <- Conditions[j]
        # print(cond2)
        compare <- paste0(cond1," vs ", cond2)
        if (! compare %in% dontcompare ){print(compare)}
        #  newcol <- ""
        #  print(newcol_reverse)
      }
    }
}



Diff_mean_fun <- function(concat_df, dontcompare = NULL) {

  if(is.null(dontcompare)) {dontcompare <- c() } 
  # return(dontcompare)

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions
  replicate_numbers <- res$replicate_numbers


  for (i in 1:length(Conditions)) {
    for (j in (1+i):length(Conditions)) {
      if (j>length(Conditions)){ break}
      cond1 <- Conditions[i]
      cond2 <- Conditions[j]

      ## making sure to not do useless comparison for the user,
      compare <- paste0(cond1," vs ", cond2)
      if (! compare %in% dontcompare) {

     # print(paste0(cond1,"_",cond2))
     # print(paste0(cond2,"_",cond1))
     # print(cond1,cond2)
      mean_value_condi_name1 <- paste0("mean_value_", cond1)
      mean_value_condi_name2 <- paste0("mean_value_", cond2)

      mean_Fx_condi_name1 <- paste0("mean_Fx_", cond1)
      mean_Fx_condi_name2 <- paste0("mean_Fx_", cond2)

      Diff_meanValue_name1 <- paste0("Diff_meanValue_",cond1,"_",cond2) 
      Diff_meanValue_name2 <- paste0("Diff_meanValue_",cond2,"_",cond1) 

      Diff_meanFx_name1 <- paste0("Diff_meanFx_",cond1,"_",cond2)
      Diff_meanFx_name2 <- paste0("Diff_meanFx_",cond2,"_",cond1) 
      #comparison_result <- my_list[[element1]] - my_list[[element2]]

      concat_df[[Diff_meanValue_name1]] <- concat_df[[mean_value_condi_name1]] - concat_df[[mean_value_condi_name2]] # nolint
      concat_df[[Diff_meanValue_name2]] <- concat_df[[mean_value_condi_name2]] - concat_df[[mean_value_condi_name1]] # nolint

      concat_df[[Diff_meanFx_name1]] <- concat_df[[mean_Fx_condi_name1]] - concat_df[[mean_Fx_condi_name2]] # nolint
      concat_df[[Diff_meanFx_name2]] <- concat_df[[mean_Fx_condi_name2]] - concat_df[[mean_Fx_condi_name1]] # nolint
      }
    }
  }
  return(concat_df=concat_df)
}

modify_p_values <- function(col) {
  col <- ifelse(col == 0.000000e+00, 10^(-16), col)
  return(col)
}

dAUC_allcondi_fun <- function(concat_df, window_number, dontcompare) {

  if(is.null(dontcompare)) {dontcompare <- c()} 

  dAUC_allcondi <- concat_df  %>% 
    filter(window==round(window_number/2))  %>%
    mutate(window_size = abs(coor2-coor1), .keep = "all") %>%
    select("transcript", "gene", "strand", "window_size") %>% distinct()

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

  for (i in 1:length(Conditions)) {
    for (j in (1+i):length(Conditions)) {
      if (j > length(Conditions)){ break}
      cond1 <- Conditions[i]
      cond2 <- Conditions[j]

      ## making sure to not do useless comparison for the user,
      compare <- paste0(cond1, " vs ", cond2)
      if (! compare %in% dontcompare) {
        mean_Fx_condi_name1 <- paste0("mean_Fx_", cond1) #e,g Control
        mean_Fx_condi_name2 <- paste0("mean_Fx_", cond2) #e.g HS

        Diff_meanFx_name1 <- paste0("Diff_meanFx_",cond1,"_",cond2)
        Diff_meanFx_name2 <- paste0("Diff_meanFx_",cond2,"_",cond1) 

        df_name <- paste0("gene_summary_AUC_", Diff_meanFx_name1)
        assign(df_name,data.frame()) 

        # Get the data frame using the dynamically generated name
        dAUC_summary_condi_df <- data.frame()
        dAUC_summary_condi_df <- get(df_name)
        dAUC <- paste0("dAUC_",Diff_meanFx_name2)
        p_dAUC <- paste0("p_dAUC_",Diff_meanFx_name2)
        D_dAUC <- paste0("D_dAUC_",Diff_meanFx_name2)

        dAUC_summary_condi_df  <- concat_df %>%
          group_by(transcript) %>% arrange(coord) %>%
          reframe(gene=gene[1], strand=strand, transcript=transcript, 
            !!dAUC := trapz(coord,!!sym(Diff_meanFx_name2)),# delta AUC
            !!p_dAUC := ks.test(!!sym(mean_Fx_condi_name1),
            !!sym(mean_Fx_condi_name2))$p.value,
            !!D_dAUC := ks.test(!!sym(mean_Fx_condi_name1),
            !!sym(mean_Fx_condi_name2))$statistic,) %>% dplyr::distinct()
        assign(df_name, dAUC_summary_condi_df)
        dAUC_allcondi <- left_join(dAUC_allcondi, dAUC_summary_condi_df, by = c("transcript", "gene","strand")) # nolint
      }
    }
  }

  dAUC_allcondi <- dAUC_allcondi %>%
  mutate(across(contains("p_dAUC"), ~ modify_p_values(.))) 

  # Get the column names containing "p_dAUC"
  p_dAUC_columns <- grep("p_dAUC", colnames(dAUC_allcondi), value = TRUE)

  # Loop through each column, calculate adjusted p-values, and create new columns
  for (col_name in p_dAUC_columns) {
    print(col_name)
    adjusted_col_name <- paste0("adjFDR_", col_name)
    dAUC_allcondi[[adjusted_col_name]] <- p.adjust(dAUC_allcondi[[col_name]], method = "fdr")
  }
  return(dAUC_allcondi)
}



# Calculate the Area Under Curve (AUC), All conditions vs y=x 
# Calculate Mean Value over the full gene body in All conditions.

AUC_allcondi_fun <- function(concat_df,window_number) {

  # Load a library using require()
  if (!require(pracma)) {
    install.packages("pracma")
    library(pracma)
  }


  AUC_allcondi <- concat_df %>% filter(window==round(window_number/2))  %>% mutate(window_size = abs(coor2-coor1), .keep = "all") %>% select("transcript", "gene", "strand", "window_size") %>% distinct() # nolint
  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions
  n <- window_number
  cumulative_density <- seq(1, n) / n

  for (cond in Conditions) {
    #already present columns
    diff_Fx_condi_name <- paste0("diff_Fx_", cond)
    mean_value_condi_name <- paste0("mean_value_", cond)
    mean_Fx_condi_name <- paste0("mean_Fx_", cond)

    df_name <- paste0("gene_summary_AUC_", cond)
    assign(df_name,data.frame())

    # Get the data frame using the dynamically generated name
    AUC_summary_condi_df <- data.frame()
    AUC_summary_condi_df <- get(df_name)

    # new column
    AUC <- paste0("AUC_", cond)
    p_AUC <- paste0("p_AUC_", cond)
    D_AUC <- paste0("D_AUC_", cond)
    MeanValueFull_condi_name <- paste0("MeanValueFull_", cond) #mean value over the full gene body

    AUC_summary_condi_df <- concat_df %>% group_by(transcript) %>%
      arrange(coord) %>%
      reframe(gene=gene[1], 
              !!AUC := trapz(coord,!!sym(diff_Fx_condi_name)) ,
              !!p_AUC := ks.test(!!sym(mean_Fx_condi_name),cumulative_density)$p.value, # nolint
              !!D_AUC := ks.test(!!sym(mean_Fx_condi_name),cumulative_density)$statistic, # nolint
              strand = strand, transcript=transcript,
              !!(MeanValueFull_condi_name) :=mean(!!sym(mean_value_condi_name))) %>% # nolint
      dplyr::distinct()
    assign(df_name, AUC_summary_condi_df)
    AUC_allcondi <- left_join(AUC_allcondi, AUC_summary_condi_df, by = c("transcript", "gene","strand")) # nolint
  }

  AUC_allcondi <- AUC_allcondi %>%
  mutate(across(contains("p_AUC"), ~ modify_p_values(.)))

  # Get the column names containing "p_dAUC"
  p_AUC_columns <- grep("p_AUC", colnames(AUC_allcondi), value = TRUE)

  # Loop through each column, calculate adjusted p-values, and create new columns
  for (col_name in p_AUC_columns) {
    print(col_name)
    adjusted_col_name <- paste0("adjFDR_", col_name)
    AUC_allcondi[[adjusted_col_name]] <- p.adjust(AUC_allcondi[[col_name]], method = "fdr") # nolint
  }
  return(AUC_allcondi)
}
# rm(AUC_allcondi)


countNA_fun <- function(main_table, extension, working_directory) {

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  col_names <- res$col_names
  score_columns <- grep("score$", col_names, value = TRUE)  

  ## Adding NA count
  Count_NA <- main_table %>% group_by(transcript) %>%
  dplyr::reframe(gene=gene[1], strand = strand[1], across(all_of(score_columns),
    ~ mean(., na.rm = TRUE), .names = "{.col}_mean"), #score_columns obtained in the first chunk creating densities # nolint
    across(contains("score"), ~ sum(is.na(.)), .names = "{.col}_NA")) %>%
    select(-contains("_mean_NA")) %>%
    set_names(map_chr(names(.), ~ ifelse(str_detect(.x, "score") &&
      str_detect(.x, "_NA"), str_replace(.x, "_score", "_count"), .x)))

  NA_plus <- Count_NA %>%
    filter(strand=="+") %>%
    select(gene, transcript, strand, contains("plus")) %>%
    select(gene, transcript, strand, contains("NA")) %>%
    set_names(map_chr(names(.), ~ ifelse(str_detect(.x, "plus") &&
      str_detect(.x, "_NA"), str_replace(.x, ".plus_", "_strand_"), .x)))
  NA_minus <- Count_NA %>%
    filter(strand=="-") %>%
    select(gene, transcript,strand, contains("minus")) %>%
    select(gene, transcript,strand, contains("NA")) %>%
    set_names(map_chr(names(.), ~ ifelse(str_detect(.x, "minus") &&
      str_detect(.x, "_NA"), str_replace(.x, ".minus_", "_strand_"), .x)))

  rm(Count_NA)

  NA_test <- bind_rows(NA_plus, NA_minus) %>% 
    select(gene, transcript, strand, matches("_count_NA")) %>%
    select(1:4)  %>%
    dplyr::rename(Count_NA = matches("count_NA")) # I drop the other NA columns because it is the same value for all the conditions (NA depends on blacklist and unmmapable region) # nolint

  return(NA_test)
}



KneeID_fun <- function(concat_df) {

  gene_summary_allcondi <- concat_df %>% select("transcript") %>% distinct()
  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

  for (cond in Conditions) {
    diff_Fx_condi_name <- paste0("diff_Fx_", cond)
    df_name <- paste0("gene_summary_AUC_", cond)
    assign(df_name,data.frame())

    # Get the data frame using the dynamically generated name
    gene_summary_condi_df <- data.frame()
    gene_summary_condi_df <- get(df_name)

    # Perform operations on the data frame
    # For example, add columns or rows to the data frame
    max_column_name <- paste0("max_", diff_Fx_condi_name)
    knee_column_name <- paste0("knee_AUC_", cond)

    gene_summary_condi_df <- concat_df %>%
      dplyr::group_by(transcript)  %>%
      dplyr::filter(!!sym(diff_Fx_condi_name) == max(!!sym(diff_Fx_condi_name)))  %>% # !!sym() inside the filter() to convert the variable name to a symbol and then access the column using the !! operator. # nolint
      slice_min(coord, n = 1)  %>% #if equality of difference within the same gene it takes the closest knee from the TSS # nolint
      select(transcript, coord, diff_Fx_condi_name)  %>%
      dplyr::rename(!!max_column_name := !!sym(diff_Fx_condi_name)) %>%
      dplyr::rename(!!knee_column_name := coord) ## the knee position is defined as the max difference between the ecdf and y=x curve # nolint

    # Assign the modified data frame back to the dynamically generated name
    assign(df_name, gene_summary_condi_df)

    gene_summary_allcondi <- left_join(gene_summary_allcondi, gene_summary_condi_df, by = "transcript") # nolint
  }
  return(gene_summary_allcondi)
}
#Working



Attenuation_fun <- function(AUC_KS_Knee_NA_DF, concat_df, pval,Replaced) {

  res <- getting_var_names(extension, file.path(working_directory, "bedgraphs"))
  Conditions <- res$Conditions

  Complete_summary <-  left_join(concat_df, AUC_KS_Knee_NA_DF,
    by = c("gene","transcript","strand"))

  for (cond in Conditions) {
    mean_value_condi_name <- paste0("mean_value_", cond)
    print(mean_value_condi_name)
    knee_column_name <- paste0("knee_AUC_", cond)
    Attenuation_cond <- paste0("Attenuation_", cond)
    UPmean_cond <- paste0("UP_mean_", cond)
    DOWNmean_cond <- paste0("DOWN_mean_", cond)
    AUC_KS_Knee_NA_DF[[Attenuation_cond]] <- NA

    result <- Complete_summary %>% group_by(transcript) %>% arrange(coord) %>%
      dplyr::reframe(transcript=transcript[1],
        !!sym(UPmean_cond) := mean((!!sym(mean_value_condi_name))[coord <= !!sym(knee_column_name)]), # nolint
        !!sym(DOWNmean_cond) := mean((!!sym(mean_value_condi_name))[coord >= !!sym(knee_column_name) & coord <= max(coord)])) %>% # nolint
        select(transcript,!!sym(UPmean_cond),!!sym(DOWNmean_cond), !!sym(DOWNmean_cond)) %>% distinct() # nolint

    AUC_KS_Knee_NA_DF <- left_join(AUC_KS_Knee_NA_DF,result, by=c("transcript"))
    AUC_KS_Knee_NA_DF[[Attenuation_cond]] <- 100 - AUC_KS_Knee_NA_DF[[DOWNmean_cond]]/AUC_KS_Knee_NA_DF[[UPmean_cond]]*100 # nolint
  }

  if (exists("Replaced") && !is.na(Replaced)) {
    if (Replaced != "NOT") {
      for (cond in Conditions) {
        p_AUC_cond <- paste0("p_AUC_", cond)
        print(p_AUC_cond)
        AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
          mutate(!!paste0("Attenuation_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, # nolint
          Replaced, .data[[paste0("Attenuation_", cond)]])) ## replacing the Attenuation by an inout value is KS test > at threshold # nolint
        AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
          mutate(!!paste0("knee_AUC_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, .data[[paste0("knee_AUC_", cond)]])) ## replacing the knee by NA is KS test > at threshold # nolint
      }
    }
  } else {
    for (cond in Conditions) {
      p_AUC_cond <- paste0("p_AUC_", cond)
      print(p_AUC_cond)
      AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
        mutate(!!paste0("Attenuation_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, # nolint
        .data[[paste0("Attenuation_", cond)]])) ## replacing the Attenuation by an input value if KS test > at threshold # nolint
      AUC_KS_Knee_NA_DF <- AUC_KS_Knee_NA_DF %>%
        mutate(!!paste0("knee_AUC_", cond) := ifelse(.data[[p_AUC_cond]] >= pval, NA, # nolint
        .data[[paste0("knee_AUC_", cond)]])) ## replacing the knee by NA if KS test > at threshold # nolint
    }
  }
  return(AUC_KS_Knee_NA_DF)
}



##################
# MAIN
##################

results_main_table <- main_table_read(name_table, extension,
  file.path(working_directory, "bedgraphs"), 0.1) # nolint
resultsECDF <- genesECDF(main_table = results_main_table[[1]], rounding,
    expressed_transcript_name_list = results_main_table[[2]], extension,
    working_dir = file.path(working_directory, "bedgraphs"))
saveRDS(resultsECDF, file = "/g/romebioinfo/tmp/explore/resultsECDF.rds")

concat_dfFX_res <- calculates_meanFx(resultsECDF,200) ## 200 is because each gene is divided in 200 windows # nolint
saveRDS(concat_dfFX_res, file = "/g/romebioinfo/tmp/explore/concat_dfFX_res.rds")

## !! Skipping this for the moment
condition_comparison(extension,file.path(working_directory, "bedgraphs")) ## Does not return anything
dontcompare_dtag <- c("CPSF3depleted_ctrl vs CPSF3wt_HS", "CPSF3depleted_HS vs CPSF3wt_ctrl") # nolint
condition_compared(extension,working_directory,) ## Does not return anything

concat_Diff_mean_res <- Diff_mean_fun(concat_dfFX_res)
saveRDS(concat_Diff_mean_res, file = "/g/romebioinfo/tmp/explore/concat_Diff_mean_res.rds")

## Time difference of 36.69392 secs
start_time <- Sys.time()
dAUC_allcondi_res <- dAUC_allcondi_fun(concat_Diff_mean_res, 200, dontcompare_dtag) # nolint
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(dAUC_allcondi_res, file = "/g/romebioinfo/tmp/explore/dAUC_allcondi_res.rds")

start_time <- Sys.time()
AUC_allcondi_res <- AUC_allcondi_fun(concat_Diff_mean_res, 200)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(AUC_allcondi_res, file = "/g/romebioinfo/tmp/explore/AUC_allcondi_res.rds")

start_time <- Sys.time()
count_NA_res <- countNA_fun(results_main_table[[1]], extension, working_directory)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(count_NA_res, file = "/g/romebioinfo/tmp/explore/count_NA_res.rds")

start_time <- Sys.time()
KneeID_res <- KneeID_fun(concat_Diff_mean_res)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(KneeID_res, file = "/g/romebioinfo/tmp/explore/KneeID_res.rds")

AUC_KS_Knee_NA.df <- left_join(AUC_allcondi_res, dAUC_allcondi_res,
  by = c("transcript", "gene", "strand", "window_size"))  %>% 
  left_join(., KneeID_res, by = c("transcript"))  %>% 
  left_join(., count_NA_res, by = c("gene", "transcript", "strand"))
AUC_KS_Knee_NA.df <- concat_Diff_mean_res %>% group_by(transcript) %>%
  summarise( chr=chr[1], coor1=min(coor1), coor2=max(coor2), strand=strand[1],
  gene=gene[1], transcript=transcript[1], size=coor2-coor1+1) %>%
  left_join(AUC_KS_Knee_NA.df, by=c("gene", "transcript", "strand"))
saveRDS(AUC_KS_Knee_NA.df, file = "/g/romebioinfo/tmp/explore/AUC_KS_Knee_NA.df.rds") # nolint

tst_df <- Attenuation_fun(AUC_KS_Knee_NA.df, concat_Diff_mean_res, 0.1, "NOT" ) #"NOT" (not replaced) or a number for attenuation (usually 0) or NA # nolint
saveRDS(tst_df, file = "/g/romebioinfo/tmp/explore/tst_df-1.rds")

mean_value_control_full <- "MeanValueFull_ctrl"
mean_value_stress <- "MeanValueFull_HS"
AUC_ctrl <- "AUC_ctrl"
AUC_stress <- "AUC_HS"
p_value_KStest <- "adjFDR_p_dAUC_Diff_meanFx_HS_ctrl"
p_value_theoritical<- "adjFDR_p_AUC_ctrl"
tst_df <- tst_df %>%
  mutate(Universe = ifelse(window_size > 50 & Count_NA < 20 &
    !!sym(mean_value_control_full) > 0.5 & !!sym(mean_value_stress) > 0.5 &
    !!sym(p_value_theoritical)> 0.1, TRUE, FALSE)) %>%
  relocate(Universe, .before = 1)
saveRDS(tst_df, file = "/g/romebioinfo/tmp/explore/tst_df-2.rds")

tst_df <- tst_df %>% mutate(
    Group = ifelse(Universe == TRUE & !!sym(AUC_stress) > 15 & -log10(!!sym(p_value_KStest)) >1.5, "Attenuated", NA), # nolint
    Group = ifelse(Universe == TRUE & !!sym(p_value_KStest)>0.2 & !!sym(AUC_ctrl) > -10 & !!sym(AUC_ctrl) < 15 , "Outgroup", Group) # nolint
  ) %>% relocate(Group, .before = 2)
saveRDS(tst_df, file = "/g/romebioinfo/tmp/explore/tst_df-3.rds")

