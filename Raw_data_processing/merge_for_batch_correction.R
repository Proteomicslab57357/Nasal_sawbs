library(dplyr)

merge_qcpass_features <- function(split_list, target_pattern = "QCPass_Correctable") {
  
  # 1. Identify which dataframes match the target pattern
  matching_names <- names(split_list)[grepl(target_pattern, names(split_list), ignore.case = TRUE)]
  
  # Safety check: Ensure we actually found matching dataframes
  if (length(matching_names) == 0) {
    stop("Error: No dataframes matched the QCPass criteria. Check your split_list.")
  }
  
  # 2. Extract only the matching dataframes into a sub-list
  qcpass_list <- split_list[matching_names]
  
  # 3. Stack them all together into a single master dataframe
  merged_df <- bind_rows(qcpass_list)
  
  # Print Summary
  print(sprintf("Step Complete: Merged %s 'QCPass' dataframes into a single dataset.", length(matching_names)))
  print(sprintf("Total features ready for batch correction: %s", nrow(merged_df)))
  
  return(merged_df)
}

# --- How to use it in your main pipeline ---
#data_for_qcmxp_neg <- merge_qcpass_features(split_data_neg)
#data_for_qcmxp_pos <- merge_qcpass_features(split_data_pos)
# 
# # You can then export this specific dataframe to a CSV for QCMXP:
#write.csv(data_for_qcmxp_neg, r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\normalization_trial\QCMXP_Input_Ready_neg.csv}", row.names = FALSE)
#write.csv(data_for_qcmxp_pos, r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\normalization_trial\QCMXP_Input_Ready_pos.csv}", row.names = FALSE)
