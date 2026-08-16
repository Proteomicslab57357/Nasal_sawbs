library(dplyr)

blank_subtraction <- function(df) {
  col_names <- colnames(df)
  
  # 1. Protect the metadata columns (First 22 + Average_Error_PPM)
  # We only want to search for blanks in the actual intensity columns
  metadata_cols <- col_names[1:22]
  if ("Average_Error_PPM" %in% col_names) {
    metadata_cols <- c(metadata_cols, "Average_Error_PPM")
  }
  
  # Define the columns that contain the actual sample intensities
  sample_cols <- setdiff(col_names, metadata_cols)
  
  # 2. Find the word "Blank" in sample names to identify the blanks group
  blank_cols <- sample_cols[grepl("Blank", sample_cols, ignore.case = TRUE)]
  
  # 3. Identify all other actual samples (HL-60, K562, QCs, etc.)
  actual_samples <- setdiff(sample_cols, blank_cols)
  
  # Safety check: If no blanks are found, skip to prevent crashing
  if (length(blank_cols) == 0) {
    print("Warning: No columns containing 'Blank' were found. Skipping subtraction.")
    return(df)
  }
  
  # Ensure all intensity columns are numeric so we can do math on them
  df <- df %>% mutate(across(all_of(c(blank_cols, actual_samples)), ~ suppressWarnings(as.numeric(.))))
  
  # 4. Calculate the average background noise (mean of the blanks) for each feature
  df$Blank_Mean <- rowMeans(df[, blank_cols, drop = FALSE], na.rm = TRUE)
  
  # If a feature is completely missing in all blanks, rowMeans returns NaN. 
  # We convert NaN to 0 so we don't accidentally subtract "NA" from the real samples.
  df$Blank_Mean[is.nan(df$Blank_Mean) | is.na(df$Blank_Mean)] <- 0
  
  # 5. Subtract the Blank Mean from every actual sample
  for (col in actual_samples) {
    df[[col]] <- df[[col]] - df$Blank_Mean
    
    # SCIENTIFIC RULE: If intensity becomes 0 or negative, it was just background noise.
    # We convert it to NA so it is treated as a missing value in later steps.
    df[[col]] <- ifelse(df[[col]] <= 0, NA, df[[col]])
  }
  
  # 6. Clean up: Remove the original blank columns and the temporary Blank_Mean column
  df <- df %>% dplyr::select(-all_of(blank_cols), -Blank_Mean)
  
  print(sprintf("Step Complete: Blank subtraction applied using %s blank samples. Blank columns removed.", length(blank_cols)))
  
  return(df)
}

# --- How to use it ---
#data_no_blanks_neg <- blank_subtraction(data_unduplicated_neg)
#data_no_blanks_pos <- blank_subtraction(data_unduplicated_pos)
