library(dplyr)
library(tidyr)
library(openxlsx)

split_features_advanced <- function(df, metadata, presence_threshold = 0.50) {
  
  # 1. Identify sample columns dynamically
  sample_cols <- colnames(df)[grepl("IDA", colnames(df)) & !grepl("ErrorPPM", colnames(df))]
  
  # 2. Convert to Long Format (Strictly using Alignment ID)
  df_long <- df %>%
    dplyr::select(`Alignment ID`, all_of(sample_cols)) %>%
    pivot_longer(cols = -`Alignment ID`, names_to = "Sample_ID", values_to = "Intensity") %>%
    left_join(metadata, by = "Sample_ID")
  
  # 3. Identify the true biological groups dynamically
  bio_groups <- unique(df_long$Group[!grepl("QC", df_long$Group, ignore.case = TRUE) & !is.na(df_long$Group)])
  n_groups <- length(bio_groups)
  
  # ---------------------------------------------------------
  # CHECK 1: BIOLOGICAL PRESENCE (>= 40% rule)
  # ---------------------------------------------------------
  bio_presence <- df_long %>%
    filter(Group %in% bio_groups) %>%
    group_by(`Alignment ID`, Group) %>%
    summarise(
      Presence_Ratio = sum(!is.na(Intensity)) / n(),
      Is_Present = Presence_Ratio >= presence_threshold,
      .groups = "drop"
    )
  
  # Create the text string of which specific groups passed
  passing_groups_df <- bio_presence %>%
    filter(Is_Present == TRUE) %>%
    group_by(`Alignment ID`) %>%
    summarise(
      Passing_Bio_Groups = paste(Group, collapse = ", "),
      .groups = "drop"
    )
  
  # Count the groups and classify
  bio_check <- bio_presence %>%
    group_by(`Alignment ID`) %>%
    summarise(
      Groups_Present = sum(Is_Present),
      .groups = "drop"
    ) %>%
    mutate(
      Bio_Category = case_when(
        Groups_Present == n_groups ~ "1_All_Groups",
        Groups_Present == 1 ~ "2_Single_Group",
        Groups_Present == n_groups - 1 ~ "3_All_Except_One",
        Groups_Present > 1 & Groups_Present < n_groups - 1 ~ "4_Shared_Multiple",
        Groups_Present == 0 ~ "6_Discarded_Low_Presence"
      )
    )
  
  # ---------------------------------------------------------
  # CHECK 2: TECHNICAL QCMXP BATCH CORRECTION READINESS
  # ---------------------------------------------------------
  qc_batches <- unique(metadata$Batch[metadata$Group == "QC"])
  n_qc_batches <- length(qc_batches)
  
  # Calculate reads per batch
  qc_batch_reads <- df_long %>%
    filter(Group == "QC") %>%
    group_by(`Alignment ID`, Batch) %>%
    summarise(QC_Valid_Reads = sum(!is.na(Intensity)), .groups = "drop")
  
  # NEW STEP: Create the formatted QC detection string (e.g., "Batch1_4, Batch2_4")
  qc_summary_df <- qc_batch_reads %>%
    mutate(Batch_String = paste0("Batch", Batch, "_", QC_Valid_Reads)) %>%
    group_by(`Alignment ID`) %>%
    summarise(
      QC_Detection_Summary = paste(Batch_String, collapse = ", "),
      .groups = "drop"
    )
  
  # Check if it passes QCMXP requirements (>= 1 read in ALL batches)
  qc_check <- qc_batch_reads %>%
    group_by(`Alignment ID`) %>%
    summarise(
      Batches_With_QC = sum(QC_Valid_Reads >= 1),
      Pass_QCMXP = (Batches_With_QC == n_qc_batches), 
      .groups = "drop"
    )
  
  # ---------------------------------------------------------
  # COMBINE AND CLASSIFY
  # ---------------------------------------------------------
  master_classification <- bio_check %>%
    left_join(qc_check, by = "Alignment ID") %>%
    mutate(
      Final_Category = case_when(
        Bio_Category == "6_Discarded_Low_Presence" ~ "6_Discarded_Low_Presence",
        Pass_QCMXP == TRUE ~ paste0(Bio_Category, "_QCPass_Correctable"),
        Pass_QCMXP == FALSE ~ paste0(Bio_Category, "_QCFail_Uncorrectable")
      )
    )
  
  # ---------------------------------------------------------
  # ATTACH NEW COLUMNS TO MAIN DATAFRAME
  # ---------------------------------------------------------
  df_annotated <- df %>%
    left_join(passing_groups_df, by = "Alignment ID") %>%
    left_join(qc_summary_df, by = "Alignment ID") %>%
    mutate(
      Passing_Bio_Groups = replace_na(Passing_Bio_Groups, "None"),
      QC_Detection_Summary = replace_na(QC_Detection_Summary, "No_QC_Reads")
    ) 
  
  # Reorganize columns so both new summary columns are placed neatly before the intensities
  metadata_cols <- setdiff(colnames(df_annotated), c(sample_cols, "Passing_Bio_Groups", "QC_Detection_Summary"))
  df_annotated <- df_annotated %>%
    dplyr::select(all_of(metadata_cols), Passing_Bio_Groups, QC_Detection_Summary, all_of(sample_cols))
  
  # ---------------------------------------------------------
  # SPLIT DATAFRAME INTO THE SPECIFIC LIST
  # ---------------------------------------------------------
  unique_categories <- unique(master_classification$Final_Category)
  
  split_list <- list()
  total_features_tracked <- 0
  
  for (cat in unique_categories) {
    target_ids <- master_classification$`Alignment ID`[master_classification$Final_Category == cat]
    split_list[[cat]] <- df_annotated %>% filter(`Alignment ID` %in% target_ids)
    total_features_tracked <- total_features_tracked + length(target_ids)
  }
  
  print(sprintf("Step Complete: Features split based on %.0f%% presence threshold across %s biological groups.", presence_threshold * 100, n_groups))
  print(sprintf("Input Features: %s | Tracked Features Output: %s", nrow(df), total_features_tracked))
  
  for (cat in names(split_list)) {
    print(sprintf(" - %s: %s features", cat, nrow(split_list[[cat]])))
  }
  
  return(split_list)
}

# --- How to use it ---

#metadata_pos <- read.xlsx(r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\input\metadatapos.xlsx}")
#metadata_neg <- read.xlsx(r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\input\metadataneg.xlsx}")
#split_data_pos <- split_features_advanced(data_no_blanks_pos, metadata_pos, presence_threshold = 0.50)
#split_data_neg <- split_features_advanced(data_no_blanks_neg, metadata_neg, presence_threshold = 0.50)

library(openxlsx)
library(dplyr)

export_to_multisheet_excel <- function(split_list, output_file_path) {
  
  # 1. Excel has a strict 31-character limit for sheet names.
  # We need to safely truncate the names in our list if they are too long.
  original_names <- names(split_list)
  safe_names <- substr(original_names, 1, 31)
  
  # Ensure truncation didn't accidentally create duplicate sheet names
  safe_names <- make.unique(safe_names, sep = "_")
  names(split_list) <- safe_names
  
  # 2. Safety check: Remove any dataframes that ended up completely empty (0 rows)
  # Excel will sometimes throw an error if you try to write a sheet with no data.
  non_empty_list <- split_list[sapply(split_list, nrow) > 0]
  
  if (length(non_empty_list) == 0) {
    print("Warning: All dataframes are empty. No Excel file was created.")
    return(invisible(NULL))
  }
  
  # 3. Write the entire list to a single Excel file
  # openxlsx automatically maps each item in the list to its own sheet
  write.xlsx(
    x = non_empty_list, 
    file = output_file_path, 
    asTable = FALSE,      # Keeps the format clean without applying Excel "Table" styles
    overwrite = TRUE      # Allows you to overwrite the file if you run the script multiple times
  )
  
  print(sprintf("Step Complete: Exported %s non-empty sheets into '%s'.", length(non_empty_list), output_file_path))
  
  # Print the exact sheet names that were created for your records
  print("Sheets created:")
  print(names(non_empty_list))
}

# --- How to use it in your main pipeline ---
#export_to_multisheet_excel(split_data_neg, r"{D:\CCHE57357\Dr. Samah NS\REsults\final_preprocessing\input\Filtered_Features_Split_50NA_neg.xlsx}")
#export_to_multisheet_excel(split_data_pos, r"{D:\CCHE57357\Dr. Samah NS\REsults\final_preprocessing\input\Filtered_Features_Split_50NA_pos.xlsx}")
