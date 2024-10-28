# getDMPs
functions I use frequently to automate DMP analysis with the limma R package

## example use (look for differentially methylated positions (DMPs) between samples with and without diabetes)
```#remove samples with no diabetes data from variable data frame
samples_with_no_diabetes_data <- variable_dataframe$sample_name[is.na(variable_dataframe$diabetes_status)]
diabetes_dataframe <- variable_dataframe[!(variable_dataframe$sample_name %in% samples_with_no_diabetes_data)]

#remove samples with no diabetes data from beta matrix
diabetes_beta_matrix <- beta_matrix[,!(colnames(beta_matrix) %in% samples_with_no_diabetes_data)]
#before running DMP analysis probes in SNPs or on the X and Y chromosome should also be removed

#import illumina manifest (optional)
manifest <- read.csv('path-to-illumina-epicv2-manifest.csv')

#run limma analysis
library(limma)
library(dplyr)
diabetes_DMPs <- get_cat_DMPs(cat_var = 'diabetes_status',
                              var_levels = c('Yes','No'),
                              s_sheet = diabetes_dataframe,
                              beta_matrix = diabetes_beta_matrix,
                              adj_var = c('age','sex',
                                          'smoking_status'), #etc...
                              annotate_with = manifest, #or if you don't have a manifest file just leave this as NULL
                              var_type = 'categorical',
                              return_bay = FALSE)

#explore results
nrow(diabetes_DMPs$Bonferroni) #this will print the number of DMPs when using more stringent Bonferroni correction
nrow(diabetes_DMPs$FDR) #this will print the number of DMPs when using less stringent false discovery rate <5% correction

#names and details of DMPs are contained in diabetes_DMPs$Bonferroni and diabetes_DMPs$FDR
#view these data frames for more info on them```
