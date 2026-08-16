library(dplyr)
library(readxl)

filter_and_summarize_ppm <- function(df, ppm_limit = 20) {
  
  col_names <- colnames(df)
  ppm_cols <- col_names[grepl("ErrorPPM", col_names, ignore.case = TRUE)]
  
  # Ensure all PPM columns are treated as numbers
  df <- df %>% mutate(across(all_of(ppm_cols), ~ suppressWarnings(as.numeric(.))))
  
  # ---------------------------------------------------------
  # 1. CELL-LEVEL FILTRATION (Evaluate each single read)
  # ---------------------------------------------------------
  for (ppm_col in ppm_cols) {
    ppm_idx <- which(colnames(df) == ppm_col)
    intensity_col <- colnames(df)[ppm_idx + 1]
    
    df[[intensity_col]] <- suppressWarnings(as.numeric(df[[intensity_col]]))
    
    # Identify individual reads that are missing or exceed the limit (+ or -)
    invalid_ppm_mask <- is.na(df[[ppm_col]]) | abs(df[[ppm_col]]) > ppm_limit
    
    # Nullify ONLY the specific intensity read that failed
    df[[intensity_col]][invalid_ppm_mask] <- NA
    
    # Nullify the bad PPM value as well so it doesn't skew our tracking metrics later
    df[[ppm_col]][invalid_ppm_mask] <- NA
  }
  
  # ---------------------------------------------------------
  # 2. CLEANUP & METADATA RECALCULATION
  # ---------------------------------------------------------
  # Recalculate the Average Error PPM strictly using the valid reads that survived
  df$Average_Error_PPM <- rowMeans(abs(df[, ppm_cols, drop = FALSE]), na.rm = TRUE)
  
  metadata_cols <- colnames(df)[1:22]
  intensity_cols <- setdiff(colnames(df), c(metadata_cols, "Average_Error_PPM", ppm_cols))
  
  # Safety Net: If a feature (row) had bad PPMs across EVERY single sample, 
  # it is now 100% NAs. We drop these entirely dead rows.
  df$Valid_Reads <- rowSums(!is.na(df[, intensity_cols]), na.rm = TRUE)
  
  df <- df %>% 
    dplyr::filter(Valid_Reads > 0) %>%
    dplyr::select(all_of(metadata_cols), Average_Error_PPM, all_of(intensity_cols))
  
  print(sprintf("Step Complete: Evaluated each single read independently. Nullified individual reads exceeding PPM limit of +/- %s.", ppm_limit))
  
  return(df)
}

# --- How to run it ---
# data_ppm_filtered <- filter_and_summarize_ppm(raw_data, ppm_limit = 20)

# --- How to use it ---
#raw_data_pos <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_pos2.xlsx}")
#raw_data_neg <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_Neg2.xlsx}")
#errorppm_filtered_data_pos <- filter_and_summarize_ppm(norm_Data_pos, ppm_limit = 20)
#errorppm_filtered_data_neg <- filter_and_summarize_ppm(norm_Data_neg, ppm_limit = 20)

