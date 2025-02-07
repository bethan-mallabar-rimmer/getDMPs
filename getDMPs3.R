#function to get run linear model and get DMPs for categorical or continuous variables
#---------------------
#by default this function assumes that categorical variables have only
#2 levels, e.g. diabetes_status = 'Yes' or 'No'
#it is possible to adapt for 3+ levels if necessary. The code for this
#is contained in the function in case you need to adapt it but I've
#commented it out because it's messy and horrible :/

getDMPs <- function(cat_var = 'diabetes_status',
                         var_levels, #leave null for continuous, replace for categorical
                         s_sheet,
                         beta_matrix,
                         adj_var, #leave null if not adjusting for anything, otherwise use e.g. adj_var = c('age','sex','smoking_status') or adj_var = 'age'
                         annotate_with = NULL,
                         var_type = 'categorical',
                         return_bay = FALSE) {
  
  #FUNCTION INPUT:
  #cat_var is the variable being analysed e.g. diabetes_status.
  #needs have the same name as whichever column in s_sheet
  #contains this data
  
  #var_levels = the 2 groups of cat_var to find DMPs between. Only required
  #for categorical analysis. E.g. var_levels = c('Yes','No') finds DMPs
  #between samples with diabetes_status = 'Yes' vs diabetes_status = 'No'
  #must be formatted the same as data in s_sheet. E.g. if s_sheet$diabetes_status
  #= c('diabetes','no diabetes','no diabetes','diabetes') then var_levels =
  #c('diabetes','no_diabetes').
  #Order does matter! If cat_var = c('Yes','No') then hypermethylated DMPs have
  #increased methylation in samples *with* diabetes. If cat_var = c('No','Yes')
  #then hypermethylated DMPs are hypermethylated in samples *without* diabetes.
  
  #s_sheet = the data frame which contains all relevant variables as
  #columns, with one row per sample. column names must include
  #both cat_var and everything in adj_var. any samples which do not
  #have variable data (e.g. diabetes_status = NA) should be removed
  #from s_sheet before running the function
  
  #beta_matrix = your beta matrix (can also be an M value matrix),
  #sample names in columns, sites in rows. Any samples which do not
  #have variable data (e.g. diabetes_status = NA) should be removed
  #from beta_matrix before running the function
  
  #adj_var = optional, a list of variables to adjust for in the model
  #e.g. adj_var = c('age','sex','smoking_status'). Any variables in this
  #list need to be in s_sheet.
  
  #annotate_with: optional! prior to running this function, you can
  #import the Illumina EPICv2 manifest as a data frame, and put
  #annotate_with = the name of this data frame. This will automatically
  #annotate all DMPs to come out of the analysis with info in the manifest
  
  #var_type: either 'categorical' or 'continuous' 
  
  #return_bay: ignore and leave as false - was just using this to
  #return results at an earlier stage of the pipeline and fix
  #some stuff.

  if (var_type == 'continuous') {var_levels = NULL}
  
  if(ncol(beta_matrix) != nrow(s_sheet)) {
  stop('beta matrix and s_sheet must contain same number of samples')}
  if(ncol(beta_matrix) != sum(!is.na(s_sheet[,colnames(s_sheet) == cat_var]))) {
    stop(paste('some samples in s_sheet do not have a value for',cat_var,'- please
               remove missing samples'))}
  
  if (var_type == 'categorical') {
    
    print('convert category variable into a factor')
    cat_f <- as.factor(s_sheet[,colnames(s_sheet) == cat_var])
    
    print('make model matrix to input into linear model')
    if (is.null(adj_var)) {
      ds <- model.matrix(as.formula('~0 + cat_f'))
    } else {
      ds <- model.matrix(as.formula(
      paste0('~0 + cat_f + s_sheet$', paste(adj_var, collapse=" + s_sheet$"))
    ))
    }
    
    #sometimes the above step removes samples, even when those samples don't
    #have any missing data - not sure why
    #at this point make sure ds and beta_matrix match by removing the same
    #samples from both and print a warning
    if (nrow(ds) != length(cat_f)) {
      x <- (1:ncol(beta_matrix))[!(1:ncol(beta_matrix) %in% rownames(ds))]
      print('warning: limma removed some samples during this stage - not sure why this sometimes happens')
      print(paste('These were',paste(colnames(beta_matrix)[x], collapse=', ')))
      print('continuing DMP analysis without these samples.')
      beta_matrix <- beta_matrix[,-x]
    }
    
    
    print('run linear model')
    tlm <- lmFit(beta_matrix, ds)
    
    print('rename columns to avoid errors')
    temp_dlm <- col_rename(column_names = adj_var,
                           design_sample = ds,
                           linear_model = tlm,
                           ss = s_sheet,
                           cat_factor = cat_f)
    tlm2 <- temp_dlm[[1]]
    ds2 <- temp_dlm[[2]]
    
    print('make contrast matrix')

    contm <- makeContrasts(case_vs_control = paste0('cat_f',var_levels[1],' - cat_f',var_levels[2]), levels=ds2)

    print('fit linear model to contrasts')
    clm <- contrasts.fit(tlm2, contm)
      
    print('empirical Bayesian estimation of differentially expressed genes (DEGs)')
    bay <- eBayes(clm)
      
    if (return_bay == TRUE) {
      return(bay)
    }
      
    print('build ranked list of DMPs for each comparison')
    topt_bon <- topTable(bay, adjust.method="bonferroni", number=nrow(beta_matrix))
    topt_fdr <- topTable(bay, adjust.method="fdr", number=nrow(beta_matrix))
    topt_bon <- topt_bon[topt_bon$adj.P.Val<0.05,]
    topt_fdr <- topt_fdr[topt_fdr$adj.P.Val<0.05,]
    
    if (!is.null(annotate_with)) {
      print('annotate with manifest')
      if (!('IlmnID' %in% colnames(annotate_with))) {
        print('Warning: IlmnID column name not in manifest! Skipping annotation step.')
      } else {
        DMPs_bon <- merge(topt_bon, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
        DMPs_fdr <- merge(topt_fdr, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      }
    } else {
      DMPs_bon <- topt_bon
      DMPs_fdr <- topt_fdr
    }
      
    DMPs <- list("Bonferroni"=DMPs_bon, "FDR"=DMPs_fdr)
    return(DMPs)
      
    #this commented out bit of code can be adapted to look for DMPs
    #in a variable with 3+ levels
    #e.g. in this case it compares tumours with fast vs intermediate vs
    #slow growth rate.
    #the only way to do this was with 3 comparisons: fast v intermediate,
    #intermediate v slow, and fast v slow.
    #so if using this code, it unfortunately needs to manually adapted
    #to the variable you are looking at
    
    #else if (cat_var == 'growth_rate_category3') {
      #contm1 <- makeContrasts(fast_vs_intermediate = 'cat_fFast - cat_fIntermediate',levels = ds2)
      #contm2 <- makeContrasts(intermediate_vs_slow = 'cat_fIntermediate - cat_fSlow',levels = ds2)
      #contm3 <- makeContrasts(fast_vs_slow = 'cat_fFast - cat_fSlow',levels = ds2)
      
      #print('fit linear model to contrasts')
      #clm1 <- contrasts.fit(tlm2, contm1)
      #clm2 <- contrasts.fit(tlm2, contm2)
      #clm3 <- contrasts.fit(tlm2, contm3)
      
      #print('empirical Bayesian estimation of differentially expressed genes (DEGs)')
      #bay1 <- eBayes(clm1)
      #bay2 <- eBayes(clm2)
      #bay3 <- eBayes(clm3)
      
      #bay <- list(fast_vs_intermediate = bay1, intermediate_vs_slow = bay2, fast_vs_slow = bay3)
      #if (return_bay == TRUE) {
      #  return(bay)
      #}
      
      #print('build ranked list of DMPs for each comparison')
      #topt_bon1 <- topTable(bay1, adjust.method="bonferroni", number=nrow(beta_matrix))
      #topt_fdr1 <- topTable(bay1, adjust.method="fdr", number=nrow(beta_matrix))
      #topt_bon1 <- topt_bon1[topt_bon1$adj.P.Val<0.05,]
      #topt_fdr1 <- topt_fdr1[topt_fdr1$adj.P.Val<0.05,]
      #topt_bon2 <- topTable(bay2, adjust.method="bonferroni", number=nrow(beta_matrix))
      #topt_fdr2 <- topTable(bay2, adjust.method="fdr", number=nrow(beta_matrix))
      #topt_bon2 <- topt_bon2[topt_bon2$adj.P.Val<0.05,]
      #topt_fdr2 <- topt_fdr2[topt_fdr2$adj.P.Val<0.05,]
      #topt_bon3 <- topTable(bay3, adjust.method="bonferroni", number=nrow(beta_matrix))
      #topt_fdr3 <- topTable(bay3, adjust.method="fdr", number=nrow(beta_matrix))
      #topt_bon3 <- topt_bon3[topt_bon3$adj.P.Val<0.05,]
      #topt_fdr3 <- topt_fdr3[topt_fdr3$adj.P.Val<0.05,]
      
      #print('annotate with manifest')
      #DMPs_bon1 <- merge(topt_bon1, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      #DMPs_fdr1 <- merge(topt_fdr1, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      #DMPs_bon2 <- merge(topt_bon2, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      #DMPs_fdr2 <- merge(topt_fdr2, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      #DMPs_bon3 <- merge(topt_bon3, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      #DMPs_fdr3 <- merge(topt_fdr3, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      
      #DMPs <- list("fast_vs_intermediate"=list("Bonferroni"=DMPs_bon1, "FDR"=DMPs_fdr1),
       #            "slow_vs_intermediate"=list("Bonferroni"=DMPs_bon2, "FDR"=DMPs_fdr2),
        #           "fast_vs_slow"=list("Bonferroni"=DMPs_bon3, "FDR"=DMPs_fdr3))
      #return(DMPs)
    #}
    
    
  } else if (var_type == 'continuous') {
    
    cat_f <- s_sheet[,colnames(s_sheet) == cat_var]
    
    print('make model matrix to input into linear model')
    if (is.null(adj_var)) {
      ds <- model.matrix(as.formula('~cat_f'))
    } else {
      ds <- model.matrix(as.formula(
      paste0('~cat_f + s_sheet$', paste(adj_var, collapse=" + s_sheet$"))
    ))
    }
    
    print('run linear model')
    tlm <- lmFit(beta_matrix, ds)
    
    print('empirical Bayesian estimation of differentially expressed genes (DEGs)')
    bay <- eBayes(tlm)
    
    if (return_bay == TRUE) {
      return(bay)
    }
    
    print('build ranked list of DMPs for each comparison')
    topt_bon <- topTable(bay, adjust.method="bonferroni", number=nrow(beta_matrix))
    topt_fdr <- topTable(bay, adjust.method="fdr", number=nrow(beta_matrix))
    topt_bon <- topt_bon[topt_bon$adj.P.Val<0.05,]
    topt_fdr <- topt_fdr[topt_fdr$adj.P.Val<0.05,]
    
    if (!is.null(annotate_with)) {
      print('annotate with manifest')
      DMPs_bon <- merge(topt_bon, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
      DMPs_fdr <- merge(topt_fdr, annotate_with, by.x="row.names", by.y="IlmnID", sort=F)
    } else {
      DMPs_bon <- topt_bon
      DMPs_fdr <- topt_fdr
    }
    
    DMPs <- list("Bonferroni"=DMPs_bon, "FDR"=DMPs_fdr)
    return(DMPs)
  }
}


#function to rename columns of design sample and linear model for Limma analysis.
#This is used by the main function - it helps prevent an error where lmFit gives
#the data column names that R is not happy with e.g. containing spaces or starting with numbers
col_rename <- function(column_names,
                       linear_model,
                       design_sample,
                       ss,
                       cat_factor) {
  nlev <- length(levels(cat_factor))
  x <- 0
  y <- c()
  for (i in 1:length(column_names)) {
    if (is.character(ss[,colnames(ss) == column_names[i]]) |
        is.factor(ss[,colnames(ss) == column_names[i]])) {
      z <- levels(as.factor(ss[,colnames(ss) == column_names[i]]))[-1]
      z <- gsub(' ','_',z)
      x <- x + length(levels(as.factor(ss[,colnames(ss) == column_names[i]]))) - 1
      y <- c(y, z)
    } else {
      x <- x + 1
      y <- c(y, column_names[i])
    }
  }
  for (i in 1:length(y)) {
    dimnames(linear_model$cov.coefficients)[[2]][[(i+nlev)]]<- y[i]
    dimnames(linear_model$coefficients)[[2]][[(i+nlev)]]<- y[i]
    dimnames(linear_model$stdev.unscaled)[[2]][[(i+nlev)]]<- y[i]
    dimnames(linear_model$design)[[2]][[(i+nlev)]]<- y[i]
    colnames(design_sample)[(i+nlev)] <- y[i]
  }
  return(list(linear_model, design_sample))
}
