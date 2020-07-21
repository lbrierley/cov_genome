#######################
# Load packages, data #
#######################

library(caret)
library(e1071)
library(janitor)
library(matrixStats)
library(magrittr)
library(pROC)
library(randomForest)
library(tidyverse)

# Load in previous ML data
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

# # Downsample groups of outcome variable
# # For each group, if n > 200 take proportional sample based on coronavirus taxids represented to give approximately 200 data points
# templist = list() # Setup empty list
# 
# for (i in 1:nlevels(model_df_predownsample$outcome)) {
#   
#   temp <- model_df_predownsample %>%
#     filter(group_name == levels(model_df_predownsample$outcome)[i]) 
#   
#   if (nrow(temp) > 200){
#     
#     # Anything less than nrow(temp)/200, take it, otherwise randomly sample
#     temp %<>% slice(suppressWarnings(createDataPartition(temp$taxid, p = 200/nrow(temp), list = FALSE, times=1)))
#     
#   }
#   
#   templist[[i]] <- temp
#   
# }
# 
# model_df <- bind_rows(templist)

# Downsample group:species of outcome variable
# For each group, if n > 20 downsample to 20
templist = list() # Setup empty list

model_df_predownsample %<>% mutate(sampler = factor(paste(outcome,childtaxa_name)))

for (i in 1:nlevels(model_df_predownsample$sampler)) {
  
  temp <- model_df_predownsample %>%
    filter(sampler == levels(model_df_predownsample$sampler)[i]) 
  
  if (nrow(temp) > 20){
    
    # Randomly sample down to 20
    temp %<>% sample_n(20)
    
  }
  
  templist[[i]] <- temp
  
}

model_df <- bind_rows(templist) %>% select(-sampler)

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
nloops <- 1
training_rows = createDataPartition(paste(model_df$outcome, model_df$genus), p = s, list = TRUE, times=nloops)
model_df %<>% select(-genus) %>% mutate(outcome = as.factor(outcome))

# Specify formula used
preds <- model_df %>% select(-outcome, -taxid, -childtaxa_name, -accessionversion) %>% remove_constant %>% names
formula_used <- as.formula(paste0("outcome ~ ", paste(preds, collapse='+')))

# Raw probability of outcome classes
vector_probs <- prop.table(table(model_df$outcome))

###########
# Run SVM
###########

full_svm_analysis <- function(x, data){
  
  train <- data %>%
    filter(row_number() %in% x) 
  
  test <- data %>%
    filter(!(row_number() %in% x))
  
  test_names <- data %>%
    filter(!(row_number() %in% x)) %>%
    select(taxid, childtaxa_name, accessionversion)
  
  
  ### TEST AREA CARET
  # Train and validate SVM (tuning gamma and cost parameters) through 10-fold cross-validation using caret (not working)
  train(x = train %>% select(preds) %>% as.data.frame,
        y = train %>% pull(outcome),
        method = "svmRadial",
        type="C-svc",
        preProc = c("center", "scale"),
        metric = "Accuracy",
        trControl = trainControl(method = 'cv',
                     number = 10,
                     verboseIter = TRUE,
                     classProbs = TRUE),
        tuneGrid = expand.grid(.sigma = seq(from = 0.01, to = 1, length = 3),
                               .C = seq(from = 0.01, to = 1, length = 3))
  )
  
  
  # Train and validate SVM (tuning gamma and cost parameters) through 10-fold cross-validation using caret (not working)
  train(x = train %>% select(preds) %>% as.matrix,
        y = train %>% pull(outcome),
        method = "rf",
        preProc = c("center", "scale"),
        metric = "Accuracy",
    #    nodesize = ,
    #    ntree = ,
        trControl = trainControl(method = 'cv',
                     number = 10,
                     verboseIter = TRUE,
                     classProbs = TRUE),
        tuneGrid = expand.grid(
        # .nodesize = seq(from = 5, to = 20, length = 5),
        #  .ntree =  seq(from = 500, to = 2000, length = 5),  
          .mtry = seq(from = 5, to = 20, length = 5))
  )
  
  ####
  
  # Train and validate SVM (tuning gamma and cost parameters) through 10-fold cross-validation using e1071
  tuned_svm <- tune.svm(formula_used, 
                        data=train,
                        type="C-classification",
                        kernel="radial",
                        gamma = seq(from = 0.01, to = 1, length = 2),
                        cost = seq(from = 0.01, to = 1, length = 2),
                        probability=TRUE,
                        tune_control = tune.control(cross = 10)) 
  
  # Obtain best performing model
  svm <- tuned_svm$best.model
  
  # Predictions for test set
  predict_class_test <- predict(svm, newdata=test, type="response")
  predict_prob_test <- predict(svm, newdata=test, probability=TRUE) %>% attributes() %>% .$probabilities %>% as.data.frame
  
  # Store individual SVM results as a list 
  list(
    tuned_obj = tuned_svm,
    best_params = tuned_svm$best.parameters,
    
    # varimp = varimp,
    
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
svm_list <- lapply(training_rows, full_svm_analysis, data=model_df)
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
  tuned_rf <- tune.randomForest(formula_used, 
                                data=train,
                                importance=TRUE, 
                                na.action=na.exclude, 
                                keep.forest=TRUE, 
                                proximity=TRUE, 
                                nodesize = seq(from = 5, to = 20, length = 5),
                                ntree =  seq(from = 500, to = 2000, length = 5),   # reduce this?
                                mtry = seq(from = 5, to = 20, length = 5))
  
  random_forest <- tuned_rf$best.model
  
  # Predictions for test set
  predict_class_test <- predict(random_forest, newdata=test, type="response")
  predict_prob_test <- predict(random_forest, newdata=test, type="prob")
  
  # Store individual RF results as a list
  list(
    
    tuned_obj = tuned_rf,
    best_params = tuned_rf$best.parameters,
    
    # Variable importance
    varimp = random_forest$importance %>%
      as.data.frame() %>%
      rownames_to_column("name") %>%
      mutate(relGini = MeanDecreaseGini/max(MeanDecreaseGini)),
    
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
rf_list <- lapply(training_rows, full_rf_analysis, data=model_df)
rf_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df,
     model_df_predownsample,
     rf_start, rf_end, svm_start, svm_end,
     rf_list, svm_list,
     file=paste0("listresults_ml_vector_", format(Sys.time(), "%d_%m_%y"), ".RData"))