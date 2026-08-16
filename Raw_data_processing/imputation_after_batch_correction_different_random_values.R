library(dplyr)
library(tidyr)

impute_post_batch_correction <- function(input_csv, output_csv = "Imputed_Final_Data.csv", impute_qcs = FALSE, max_na_ratio = 0.60) {
  
  # 1. Read the raw data cleanly
  raw_df <- read.csv(input_csv, header = FALSE, stringsAsFactors = FALSE)
  
  sample_names <- as.character(raw_df[1, -1])
  group_labels <- as.character(raw_df[2, -1])
  
  feature_ids <- raw_df[3:nrow(raw_df), 1]
  intensity_matrix <- raw_df[3:nrow(raw_df), -1]
  
  intensity_matrix <- as.data.frame(lapply(intensity_matrix, as.numeric))
  colnames(intensity_matrix) <- sample_names
  intensity_matrix$Feature_ID <- feature_ids
  
  # 2. Pivot to Long format for grouped math
  df_long <- intensity_matrix %>%
    pivot_longer(cols = -Feature_ID, names_to = "Sample", values_to = "Intensity")
  
  sample_to_group <- data.frame(Sample = sample_names, Group = group_labels, stringsAsFactors = FALSE)
  df_long <- df_long %>% left_join(sample_to_group, by = "Sample")
  
  print(sprintf("Evaluating NA thresholds (<= %s%%) and applying exact runif(1) imputation per group...", max_na_ratio * 100))
  
  # 3. Apply the strictly isolated Imputation Logic
  df_imputed_long <- df_long %>%
    group_by(Feature_ID, Group) %>%
    mutate(
      # Calculate median and real-time missing value percentage
      grp_med = median(Intensity, na.rm = TRUE),
      na_ratio = sum(is.na(Intensity)) / n(),
      
      # Does this specific biological group have <= 60% NAs?
      passes_cutoff = na_ratio <= max_na_ratio,
      
      # Generate the single matching random number for the group
      single_impute_val = ifelse(
        passes_cutoff & !is.na(grp_med) & grp_med > 0, 
        runif(n(), grp_med - (grp_med * 0.01), grp_med + (grp_med * 0.01)), 
        NA
      ),
      
      # Prevent QC imputation unless explicitly requested by the user
      is_valid_target = if(impute_qcs) TRUE else !grepl("QC", Group, ignore.case = TRUE),
      
      # Replace all NAs in this group with that single identical value
      Intensity = ifelse(
        is.na(Intensity) & is_valid_target & !is.na(single_impute_val),
        single_impute_val,
        Intensity
      )
    ) %>%
    ungroup() %>%
    dplyr::select(Feature_ID, Sample, Intensity)
  
  # 4. Pivot back to Wide format
  df_wide <- df_imputed_long %>%
    pivot_wider(names_from = "Sample", values_from = "Intensity")
  
  df_wide <- df_wide %>% dplyr::select(Feature_ID, all_of(sample_names))
  
  # 5. Reconstruct the original CSV structure
  row1 <- c("Sample", sample_names)
  row2 <- c("Label", group_labels)
  
  df_wide_char <- as.data.frame(lapply(df_wide, as.character), stringsAsFactors = FALSE)
  final_df <- rbind(row1, row2, df_wide_char)
  
  # 6. Save out the completed dataset
  write.table(final_df, output_csv, sep = ",", row.names = FALSE, col.names = FALSE, na = "")
  
  print(sprintf("Step Complete: Missing values imputed. Final matrix saved to '%s'.", output_csv))
  
  return(final_df)
}


#impute_post_batch_correction(input_csv = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\output\EF_DM_BS_NA\between_batch_correction_simp.csv}",
#                             output_csv = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\output\EF_DM_BS_NA\Final_Corrected_and_Imputed_different_random_values.csv}",
#                             impute_qcs = TRUE, max_na_ratio = 0.60)
