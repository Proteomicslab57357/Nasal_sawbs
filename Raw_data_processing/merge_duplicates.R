library(dplyr)

merge_duplicates_by_inchikey <- function(df) {
  
  # 1. Custom function to safely aggregate rows
  aggregate_values <- function(x) {
    if (all(is.na(x))) return(NA)
    
    # If the column is numeric (like intensities or MZ/RT), take the maximum value
    if (is.numeric(x)) {
      return(max(x, na.rm = TRUE))
    } else {
      # If it's a character column (metadata), just take the first non-null/non-NA value
      # Using max() on text is slow and can scramble metadata alphabetically
      x_valid <- x[!is.na(x) & x != "" & tolower(x) != "null"]
      if (length(x_valid) == 0) return(x[1])
      return(x_valid[1])
    }
  }
  
  # 2. Store the original column order so we can restore it exactly at the end
  original_cols <- colnames(df)
  
  # 3. Create a safe merging identifier
  df <- df %>%
    mutate(
      Merge_ID = case_when(
        # If INCHIKEY is missing or "null", use Alignment ID so unknowns aren't merged together
        is.na(INCHIKEY) | trimws(INCHIKEY) == "" | tolower(INCHIKEY) == "null" ~ as.character(`Alignment ID`),
        # Otherwise, use the valid INCHIKEY
        TRUE ~ as.character(INCHIKEY)
      )
    )
  
  # 4. Group by the safe Merge_ID and condense the duplicates
  df_merged <- df %>%
    group_by(Merge_ID) %>%
    summarise(across(everything(), aggregate_values), .groups = "drop")
  
  # 5. Remove the temporary Merge_ID column and restore original column order
  df_merged <- df_merged %>% 
    dplyr::select(all_of(original_cols))
  
  # Calculate how many duplicates were squashed
  rows_removed <- nrow(df) - nrow(df_merged)
  print(sprintf("Step Complete: Duplicates merged based on INCHIKEY. %s duplicate rows were combined.", rows_removed))
  
  return(df_merged)
}

# --- How to use it in your main pipeline ---
#data_unduplicated_neg <- merge_duplicates_by_inchikey(errorppm_filtered_data_neg)
#data_unduplicated_pos <- merge_duplicates_by_inchikey(errorppm_filtered_data_pos)
