library(openxlsx)
library(readxl)
library(dplyr)

labPQN <- function(x, qc) {
  as.numeric(x) / median(as.numeric(x) / median(as.numeric(qc), na.rm = T), na.rm = T)
}

lab_normalization <- function(data) {
  
# remove the two columns of the data containing filetype and injection order  
data <- data[-(1:2),]
# Create sample names

sampleColumns = data %>%   
  #dplyr::select(-contains("QC"), -contains("Blank")) %>%
  dplyr::select(contains("IDA")) %>%
  colnames()

sampleIndices = which(names(data) %in% sampleColumns)
errorppmIndices = sampleIndices - 1
errorPPMColumns =  names(data[, errorppmIndices])
fiTable = sort(c(errorppmIndices,sampleIndices))

# Confirming that the data is numeric
data[, sampleColumns] = apply(data[, sampleColumns], 2, as.numeric)

nsData = data %>% 
  mutate(across(all_of(sampleColumns), ~ ifelse(. == 0, NA, .))) %>% as_tibble()

# Discovering the reference sample
refSample = nsData %>% 
  dplyr::select(all_of(sampleColumns)) %>% 
  colSums(., na.rm=T) %>% 
  which.max() %>% 
  names() %>% 
  nsData[[.]]

ref_sample_name <- nsData %>% 
  dplyr::select(all_of(sampleColumns)) %>% 
  colSums(., na.rm=T) %>% 
  which.max() %>% 
  names()

print(ref_sample_name)

data_normalized = nsData %>% mutate(across(all_of(sampleColumns), ~ labPQN(., refSample)))
print('Normalization is finished')

return(data_normalized)
}

#raw_data_pos <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_pos2.xlsx}")
#raw_data_neg <- read_excel(r"{D:/CCHE57357/Dr. Samah NS/REsults/final_preprocessing/input/GNP_TABLE_Neg2.xlsx}")

#norm_Data_neg <- lab_normalization(raw_data_neg)
#norm_Data_pos <- lab_normalization(raw_data_pos)
#write.csv(norm_Data_neg, file = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\normalization_trial\normalized_data_neg.csv}")
#write.csv(norm_Data_pos, file = r"{D:\CCHE57357\Dr. Samah NS\REsults\new_trial\separate_functions\normalization_trial\normalized_data_pos.csv}")

