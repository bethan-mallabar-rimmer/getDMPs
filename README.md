# getDMPs
functions I use frequently to automate DMP analysis with the limma R package

## example use of getDMPs (look for differentially methylated positions (DMPs) between samples with and without diabetes)
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
diabetes_DMPs <- getDMPs(cat_var = 'diabetes_status',
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
#view these data frames for more info on them
```

## function inputs for getDMPs
The function can be used for categorical or continuous variables. By default the function assumes that categorical variables have only 2 levels, e.g. diabetes_status = 'Yes' or 'No'.
It is possible to adapt for 3+ levels if necessary. The code for this is contained in the function in case you need to adapt it but I've commented it out because it's messy and horrible :/

<b>`cat_var`</b> is the variable being analysed e.g. diabetes_status. Needs have the same name as whichever column in s_sheet contains this data.

<b>`var_levels`</b> = the 2 groups of cat_var to find DMPs between. Only required for categorical analysis. E.g. var_levels = c('Yes','No') finds DMPs between samples with diabetes_status = 'Yes' vs diabetes_status = 'No'. Must be formatted the same as data in s_sheet. E.g. if s_sheet$diabetes_status = c('diabetes','no diabetes','no diabetes','diabetes') then var_levels = c('diabetes','no_diabetes'). </br>
Note: <b>order does matter!</b> If cat_var = c('Yes','No') then hypermethylated DMPs have increased methylation in samples *with* diabetes. If cat_var = c('No','Yes') then hypermethylated DMPs are hypermethylated in samples *without* diabetes.
  
<b>`s_sheet`</b> = the data frame which contains all relevant variables as columns, with one row per sample. Column names must include both the argument/input to `cat_var` (e.g. diabetes_status as a column name) and everything in `adj_var`. Any samples which do not have variable data (e.g. diabetes_status = NA) should be removed
from s_sheet before running the function.
  
<b>`beta_matrix`</b> = your beta matrix (can also be an M value matrix - in fact it probably should be as this is more statistically valid), with sample names in columns, sites in rows. Any samples which do not have variable data (e.g. diabetes_status = NA) should be removed from beta_matrix before running the function
  
<b>`adj_var`</b> = a list of variables to adjust for in the model e.g. c('age','sex','smoking_status'). Any variables in this list need to be in  the column names of s_sheet
  
<b>`annotate_with`</b>: optional! prior to running this function, you can import the Illumina manifest for your platform (e.g. EPIC or EPICv2) as a data frame, and put annotate_with = the name of this data frame. This will automatically annotate all DMPs to come out of the analysis with info in the manifest. The manifest must have the column IlmnID for this to work.

<b>`var_type`</b>: either 'categorical' or 'continuous' 
  
<b>`return_bay`</b>: ignore and leave as false - I was just using this to return results at an earlier stage of the pipeline and fix some stuff.


