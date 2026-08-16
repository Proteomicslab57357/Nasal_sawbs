
# 1. Load Functions
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/labPQN_normalization.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/errorppm_filtration.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/merge_duplicates.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/Blank_subtraction.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/data_splitting.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/merge_for_batch_correction.R")
source("D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/imputation_after_batch_correction_different_random_values.R")
#source("07_export_multisheet.R")

# 2. Load Raw Data & Metadata
raw_data_neg <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_neg2.xlsx}")
metadata_neg <- read.xlsx(r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\input\metadataneg.xlsx}")
raw_data_pos <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_pos2.xlsx}")
metadata_pos <- read.xlsx(r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\input\metadatapos.xlsx}")
# 3. Run Pipeline
data_normalized_pos <- lab_normalization(raw_data_pos)
data_normalized_neg <- lab_normalization(raw_data_neg)
data_ppm_filtered_pos <- filter_and_summarize_ppm(data_normalized_pos, ppm_limit = 20)
data_ppm_filtered_neg <- filter_and_summarize_ppm(data_normalized_neg, ppm_limit = 20)
data_unduplicated_pos <- merge_duplicates_by_inchikey(data_ppm_filtered_pos)
data_unduplicated_neg <- merge_duplicates_by_inchikey(data_ppm_filtered_neg)
data_no_blanks_pos    <- blank_subtraction(data_unduplicated_pos)
data_no_blanks_neg    <- blank_subtraction(data_unduplicated_neg)
split_data_pos        <- split_features_advanced(data_no_blanks_pos, metadata_pos, presence_threshold = 0.40)
split_data_neg        <- split_features_advanced(data_no_blanks_neg, metadata_neg, presence_threshold = 0.40)

qc_correcatable_features_pos <- merge_qcpass_features(split_data_pos)
qc_correcatable_features_neg <- merge_qcpass_features(split_data_neg)

write.xlsx(qc_correcatable_features_pos, file = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\data_pos40.xlsx}")
write.xlsx(qc_correcatable_features_neg, file = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\data_neg40.xlsx}")
########################################################
# go to QCMXP to do batch correction#
########################################################
# import machine drift corrected data obtained from QCMXP
corrected_Data_pos <- "D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/data_pos_normalized_ERRFIL20_DM_BS_NA60_BC.csv"
corrected_Data_neg <- "D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/data_neg_normalized_ERRFIL20_DM_BS_NA60_BC.csv"

# output file 

imputed_data_pos <- "D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/data_pos_normalized_ERRFIL20_DM_BS_NA60_BC_imputed.csv"
imputed_data_neg <- "D:/CCHE57357/Dr. Samah NS/REsults/new_trial/separate_functions/data_neg_normalized_ERRFIL20_DM_BS_NA60_BC_imputed.csv"

#impute and output the file
impute_post_batch_correction(corrected_Data_neg, imputed_data_neg, max_na_ratio = 0.6)
impute_post_batch_correction(corrected_Data_pos, imputed_data_pos, max_na_ratio = 0.6)
