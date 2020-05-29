#######################
# Load packages, data #
#######################

library(caret)
library(e1071)
library(janitor)
library(matrixStats)
library(magrittr)
library(pROC)
library(pbapply)
library(randomForest)
library(tidyverse)

# Load in previous ML results
load("cov_ML_dfs_28_05_20.RData")

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model #
##################################################################################

# Set options
set.seed(1315)
outcome_name <- "group_name"
use_stop_codons <- TRUE
min_n_seq_category <- 50

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
if (!(outcome_name %in% names(cov_spikes_df))){
  model_df_predownsample <- cov_spikes_df %>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome_name)),
                                          by = c("taxid" = "childtaxa_id")) %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
} else {
  model_df_predownsample <- cov_spikes_df %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
}

# Downsample groups of outcome variable

# Tables of outcome variable counts
model_df_predownsample %>% 
  group_by(outcome) %>% 
  summarise(n_sequences = n(), n_species = n_distinct(taxid)) %>% 
  arrange(-n_sequences)

model_df_predownsample %>%
  filter(outcome == "aves") %>%
  filter(!is.na(genus_name)) %>%
  group_by(genus_name) %>%
  summarise(n_sequences = n(), n_species = n_distinct(taxid)) %>%
  arrange(-n_sequences)


# For each group, if n > 200 take proportional sample based on coronavirus taxids represented to give approximately 200 data points
templist = list() # Setup empty list

for (i in 1:nlevels(model_df_predownsample$outcome)) {
  
  temp <- model_df_predownsample %>%
    filter(group_name == levels(model_df_predownsample$outcome)[i]) 
  
  if (nrow(temp) > 200){
    
    # Anything less than nrow(temp)/200, take it, otherwise randomly sample
    temp %<>% slice(suppressWarnings(createDataPartition(temp$taxid, p = 200/nrow(temp), list = FALSE, times=1)))
    
  }
  
  templist[[i]] <- temp
  
}

model_df <- bind_rows(templist)

# Include nucleotide, dincucleotide and codon usage composition bias measurements, with option to exclude stop codon RSCU
if (use_stop_codons == TRUE){
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion)
} else {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias))
}

# Filter outcome to non-missing data and classes with minimum number of sequence observations
model_df %<>% group_by(outcome) %>% filter(n() >= min_n_seq_category) %>% ungroup %>% mutate(outcome = factor(outcome))

# Proportion in the training set
s <- .75 

# Set training and test sets
nloops <- 2
training_rows = createDataPartition(paste(model_df$outcome, model_df$genus), p = s, list = TRUE, times=nloops)
model_df %<>% select(-genus) %>% mutate(outcome = as.factor(outcome))

# Specify formula used
formula_used <- formula(outcome ~ .)

# Raw probability of outcome classes
prop.table(table(model_df$outcome))

###########
# Run SVM
###########

full_svm_analysis <- function(x, data){
  
  train <- data %>%
    filter(row_number() %in% x) %>%
    remove_constant %>%
    select(-c(taxid, childtaxa_name, accessionversion))
  
  test <- data %>%
    filter(!(row_number() %in% x)) %>%
    remove_constant %>%
    select(-c(taxid, childtaxa_name, accessionversion))
  
  test_names <- data %>%
    filter(!(row_number() %in% x)) %>%
    select(taxid, childtaxa_name, accessionversion)
  
  # Build SVM
  svm <- svm(formula_used, 
             data=train,
             type="C-classification",
             kernel="radial",
             cost = 1,
             probability=TRUE)  
  
  # Predictions for test set
  predict_class_test <- predict(svm, newdata=test, type="response")
  predict_prob_test <- predict(svm, newdata=test, probability=TRUE) %>% attributes() %>% .$probabilities %>% as.data.frame
  
  # Calculate relative importance of predictors - unclear how to do this for SVM
  
  # # Extract partial dependence values - not including for now as takes too long to extract
  # set <- train %>% select(-outcome) %>% names # Get predictor names
  # list_PD <- list_prob <<- vector("list", length(set)) %>% setNames(set)
  # 
  # for (i in 1:length(set)){
  #   list_PD[[i]] <- pdp::partial(svm, pred.var=set[i], train = train)
  #   list_prob[[i]] <- pdp:: partial(svm, pred.var=set[i], train = train, prob=TRUE)
  # }
  
  # Store individual SVM results as a list 
  list(
    
    # varimp = varimp,
    
    # # Partial dependence
    # list_PD = list_PD,
    # list_prob, = list_prob,
    
    # Confusion matrix, test set predictions
    matrix_test = confusionMatrix(predict_class_test, test$outcome),
    predict_test = data.frame(test_names,
                              pred_class=predict_class_test,
                              pred_prob=predict_prob_test,
                              outcome=test$outcome),
    
    
    # AUC
    AUC = multiclass.roc(response = test$outcome, predictor = predict_prob_test) %>% .$auc %>% as.numeric,
    ROC = multiclass.roc(response = test$outcome, predictor = predict_prob_test)
  )
}

svm_start <- Sys.time()
# Store ensemble results as list of lists
svm_list <- pblapply(training_rows, full_svm_analysis, data=model_df)
svm_end <- Sys.time()


######################
# Run random forests #
######################

full_rf_analysis <- function(x, data){
  
  train <- data %>%
    filter(row_number() %in% x) %>%
    remove_constant %>%
    select(-c(taxid, childtaxa_name, accessionversion))
  
  test <- data %>%
    filter(!(row_number() %in% x)) %>%
    remove_constant %>%
    select(-c(taxid, childtaxa_name, accessionversion))
  
  test_names <- data %>%
    filter(!(row_number() %in% x)) %>%
    select(taxid, childtaxa_name, accessionversion)
  
  # Build random forests
  random_forest <- randomForest(formula_used, 
                                dat=train,
                                importance=TRUE, 
                                na.action=na.exclude, 
                                keep.forest=TRUE, 
                                ntree = 5000,   # reduce this?
                                proximity=TRUE, 
                                nodesize = 10)
  
  # Predictions for test set
  predict_class_test <- predict(random_forest, newdata=test, type="response")
  predict_prob_test <- predict(random_forest, newdata=test, type="prob")
  
  
  # # Extract partial dependence values - not including for now as takes too long to extract
  # set <- train %>% select(-outcome) %>% names # Get predictor names
  # list_PD <- list_prob <<- vector("list", length(set)) %>% setNames(set)
  # 
  # for (i in 1:length(set)){
  #   list_PD[[i]] <- pdp::partial(random_forest, pred.var=set[i], train = train)
  #   list_prob[[i]] <- pdp:: partial(random_forest, pred.var=set[i], train = train, prob=TRUE)
  # }
  
  # Store individual RF results as a list
  list(
    
    # Variable importance
    varimp = random_forest$importance %>%
      as.data.frame() %>%
      rownames_to_column("name") %>%
      mutate(relGini = MeanDecreaseGini/max(MeanDecreaseGini)),
    
    # # Partial dependence
    # list_PD = list_PD,
    # list_prob = list_prob,
    
    # Confusion matrix, test set predictions
    matrix_test = confusionMatrix(predict_class_test, test$outcome),
    predict_test = data.frame(test_names,
                              pred_class=predict_class_test,
                              pred_prob=predict_prob_test,
                              outcome=test$outcome),
    
    # AUC, multiclass ROC
    AUC = multiclass.roc(response = test$outcome, predictor = predict_prob_test) %>% .$auc %>% as.numeric,
    ROC = multiclass.roc(response = test$outcome, predictor = predict_prob_test)
  )
}

rf_start <- Sys.time()
# Store ensemble RF results as list of lists
rf_list <- pblapply(training_rows, full_rf_analysis, data=model_df)
rf_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(rf_start, rf_end, svm_start, svm_end,
     rf_list, svm_list,
     file="listresults_rf_svm_12_5_20.RData")
