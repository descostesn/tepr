##works for different datasets and plot all replicates

library(tidyverse)

getting_var_names <- function(extension, working_directory) {
  # This function uses the extension and working directory to get the condition names, the number of replicates, and the variable names.
  # It needs the file names to be written in the form:
  # Condition_rep#.strand.extension such as :
  # HS_rep1.reverse.bg
  # In input: extension such as "*.bg" and the working directory
  bedgraph_files <- list.files(working_directory, pattern = extension, full.names = TRUE)
  files <- bedgraph_files %>%
    map(~{
      filename <- tools::file_path_sans_ext(basename(.))
    })
  string <- files
  var_names <- string
  
  for (i in seq_along(var_names)) {
    if (grepl("(reverse|forward)", var_names[i])) {
      var_names[i] <- gsub("reverse", "minus", var_names[i])
      var_names[i] <- gsub("forward", "plus", var_names[i])
    }
  }
  
  # Extract conditions
  Conditions <- unique(sub("(\\w+)_rep\\d+.*", "\\1", var_names)) ## verify it can work with several "_" 
  
  # Conditions <- unique(sub(".*_(\\w+)_rep\\d+.*", "\\1", var_names)) ## verify it can work with several "_" 
  # Extract replication numbers
  replicate_numbers <- unique(sub(".*_rep(\\d+).*", "\\1", var_names))

  # Setting fixed column names
  fixed_names <- c("biotype","chr", "coor1", "coor2","transcript", "gene", "strand","window","id") #biotype is added in the .tsv file
  
  num_fixed_cols <- length(fixed_names)
  # Number of variable columns based on input file
  num_var_cols <- ncol(table) - num_fixed_cols
  #Generate variable column names for category 2, alternating with category 1
  test <- rep(var_names, each=2)
  suffix <- rep(c("", "_score"), length(var_names))
  test <- paste0(test, suffix)
  col_names <- c(fixed_names, test)
  return(list(col_names=col_names,var_names=var_names, replicate_numbers=replicate_numbers, Conditions=Conditions)) 
}


# Define a function to generate color scale
generate_color_scale <- function(data, Conditions, pattern_colors, default_color = "black") {
  colors <- rep(default_color, length(data))
  # Assign colors based on patterns
  for (i in seq_along(Conditions)) {
    colors <- ifelse(grepl(Conditions[i], data), pattern_colors[i], colors)
  }
  return(colors)
}

# Ask the user for colors for each condition
ask_for_colors <- function( extension, working_directory) {
  res <- getting_var_names(extension, working_directory)
  Conditions <- res$Conditions
  replicate_numbers <- res$replicate_numbers
  pattern_colors <- vector("character", length = length(Conditions))
  for (i in seq_along(Conditions)) {
    pattern_colors[i] <- readline(prompt = paste("Enter color for pattern", Conditions[i], ": "))
  }
  return(pattern_colors)
}



#pattern_colors <- ask_for_colors(extension, working_directory)

# Define a function to plot ECDF data for a given gene
plot_gene_ECDF <- function(gene_name, resultsECDF, tst_df, extension, working_directory, pattern_colors, decimal_places) {
  res <- getting_var_names(extension, working_directory)
  Conditions <- res$Conditions
  replicate_numbers <- res$replicate_numbers
  gene_stat <- tst_df %>% filter(gene==gene_name)
  
  # Filter columns to pivot longer based on pattern
  columns_to_pivot <- paste0("Fx_", Conditions, "_rep", replicate_numbers, "_score")
  
  window_size_factor <- resultsECDF %>%
    filter(coord == 100, gene == gene_name) %>%
    mutate(window_size_factor = (coor2 - coor1) / 1000) %>%
    pull(window_size_factor)

  # Filter columns to pivot longer based on pattern for both Fx and value
  column_vector_Fx <- character()
  column_vector_value <- character()
  # Initialize a data frame to store the positions for vertical lines
  vline_data <- data.frame(Condition = character(), Knee_Value = numeric(), stringsAsFactors = FALSE)
  subtitle_text <- ""
  
  for (cond in Conditions) {
    
    # Generate names and values for AUC, KS, and knee conditions
    AUC_cond_name <- paste0("AUC_", cond)
    AUC_cond <- round(gene_stat %>% pull(AUC_cond_name), decimal_places)  # Round AUC value
    KS_cond_name <-  paste0("adjFDR_p_AUC_", cond)
    KS_cond <- formatC(gene_stat %>% pull(KS_cond_name), format = "e", digits = decimal_places)
    knee_cond_name <- paste0("knee_AUC_", cond)
    knee_cond <- round(gene_stat %>% mutate(knee_kb=!!sym(knee_cond_name)*window_size/1000) %>% pull(knee_kb))  # Round Knee value
    

    # If KS test value is less than 0.01, store the knee_cond value for the vertical line
    
    if (as.numeric(KS_cond) < 0.01) {
      vline_data <- rbind(vline_data, data.frame(Condition = cond, Knee_Value = knee_cond))
      
    } else {
      knee_cond <- NA
    }
    
    # Create a string for the current condition and append it to the subtitle
    condition_text <- paste0(cond, ": AUC = ", AUC_cond, ", KS = ", KS_cond, ", Knee (kb) = ", knee_cond)
    subtitle_text <- paste(subtitle_text, condition_text, sep = "\n")

  # list_knee <- list(AUC_cond_name, AUC_cond, KS_cond_name, KS_cond, knee_cond_name, knee_cond)
  #  print(list_knee)
   
    
    for (rep_num in replicate_numbers) {
      new_column_Fx <- paste0("Fx_", cond, "_rep", rep_num, "_score")  
      new_column_value<- paste0("value_", cond, "_rep", rep_num, "_score")
      column_vector_Fx <- c(column_vector_Fx, new_column_Fx)
      column_vector_value <- c(column_vector_value, new_column_value)
    }
  }
  
  # Pivot for Fx values
  df_long_Fx <- resultsECDF %>%
    filter(gene == gene_name) %>%
    pivot_longer(cols = all_of(column_vector_Fx),
                 names_to = "Conditions", values_to = "Fx") %>%
    mutate(Conditions = gsub("Fx_|_score", "", Conditions))
  
  # Pivot for value values
  df_long_value <- resultsECDF %>%
    filter(gene == gene_name) %>%
    pivot_longer(cols = all_of(column_vector_value),
                 names_to = "Conditions", values_to = "value") %>%
    mutate(Conditions = gsub("value_|_score", "", Conditions))
  
  # Find the common columns
  common_columns <- intersect(names(df_long_Fx), names(df_long_value))
  
  # Merge the data frames excluding Fx and value
  df_long_ECDF <- merge(df_long_Fx, df_long_value, by = common_columns)
  
  # Generate color scale
  color_scale <- generate_color_scale(df_long_ECDF$Conditions, Conditions, pattern_colors)
  
  #scaling the second y axis to prevent it to take the whole plot: 
  scale = 2*max(df_long_ECDF$value)
  
  # Plotting
  ggplot(df_long_ECDF, aes(x = coord, y = Fx, color = Conditions)) +
    geom_line(aes(x = coord * window_size_factor, y = coord / max(coord)), linetype = "dashed", color = "red") +
    geom_area(aes(x = coord * window_size_factor, y = value / scale, fill = Conditions), alpha = 0.1,linewidth=0.2, position = 'identity') +
    scale_fill_manual(values = color_scale) +
    geom_line(linewidth = 1, aes(x = coord * window_size_factor)) +
    scale_color_manual(values = color_scale) +
    
    scale_y_continuous(sec.axis = sec_axis(~ . * scale, name = "Transcription level")) +
    geom_vline(data = vline_data, aes(xintercept = Knee_Value), linetype = "dashed", color = "darkgrey") + 
    labs(x = "Distance from TSS (kb)", y = "Cumulative transcription density", title = gene_name, subtitle = subtitle_text) +
    theme_classic()
}

# Usage
# pattern_colors <- ask_for_colors( extension, working_directory )
pattern_colors <- c("#90AFBB","#FF9A04","#10AFBB", "#FC4E07")
gene_name <- "EGFR"  # Example gene
 
 my_plot <- 
  plot_gene_ECDF(gene_name, concat_dfFX_res, tst_df, extension, working_directory, pattern_colors, 2) ##last parameter is rounding

file_path <- paste("./", gene_name, ".png", sep = "")
ggsave(file_path, plot = my_plot, width = 10.5, height = 6, units = "cm", bg='transparent')


#("#00AFBB","#10AFBB","#FF9A04", "#FD1E00", "#000E00")
pattern_colors <- c("#90AFBB","#FF9A04","#10AFBB", "#FC4E07")
pattern_colors <- c("#10AFBB", "#FC4E07")

