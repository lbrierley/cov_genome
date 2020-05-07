
library(caret)
library(e1071)
library(janitor)
library(ROCR)

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model #
##################################################################################

# Prepare df for modelling
model_df <- cov_spikes_df %>% left_join(allcov_df %>% select(childtaxa_id, h_human),
                                        by = c("taxid" = "childtaxa_id")) %>% mutate(response = h_human)

# MERGE OUTCOME HELPER
if (!(outcome %in% names(df))){
  df %<>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome)),
                    by = c("taxid" = "childtaxa_id"))
}

# SELECTORS FOR EACH TYPE OF BIAS HELPER
df %>% select(matches("^[A|C|G|T][A|C|G|T]_Bias$"))
df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"))
df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias))
df %>% select(matches("^.*_aa_Bias$"))

# Proportion in the training set
s <- .75 

# Set training and test sets
set.seed(1547)
nloops <- 1
training_rows = createDataPartition(paste(cov_spikes_df$h_human, cov_spikes_df$genus), p = s, list = TRUE, times=nloops)

# Specify formula used
formula_used <- formula(response ~ .)

###########
# Run SVM
###########

full_svm_analysis <- function(x, data){
  
  train <- data %>%
    filter(row_number() %in% x) %>%
    select(response, matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>%
    remove_constant
  
  test <- data %>%
    filter(!(row_number() %in% x)) %>%
    select(response, dplyr::matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>%
    remove_constant
  
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
  predict_prob_test <- predict(svm, newdata=test, probability=TRUE) %>% attributes() %>% .$probabilities %>% as.data.frame %>% pull("1")
  predict_test_df <- data.frame(test_names,
                                pred_class=predict_class_test %>% as.character %>% as.numeric,
                                pred_prob=predict_prob_test,
                                response=test$response %>% as.character %>% as.numeric)
  
  #### ## USE THIS FOR VARIMP/CLASSIFICATION https://bgreenwell.github.io/pdp/articles/pdp-classification.html
  
  # Calculate relative importance of predictors (possible to do so with permutation algorithm, analogous to that of RF varimp)
  # varimp <- summary(gbm_ada, n.trees = opt_it) %>%  mutate(rel.inf = rel.inf/max(rel.inf)) %>% rename(varimp = rel.inf, name = var)
  
  # # Extract partial dependence values
  # list_PD <- vector("list", npd) %>% setNames(set[2:(npd+1)])
  # 
  # for (i in 2:(npd+1)){
  #   list_PD[[i-1]] <- partial(gbm_ada, pred.var=set[i], train = train, which.class="1", n.trees= opt_it)
  # }
  
  # Store individual SVM results as a list 
  list(
    
    # varimp = varimp,
    
    # # Partial dependence
    # list_PD = list_PD,
    
    # Confusion matrix, test set predictions
    matrix_test = confusionMatrix(predict_class_test, test$response, positive="1"),
    predict_test_df = predict_test_df,
    
    # AUC
    AUC = prediction(predict_test_df %>% pull(pred_class) %>% as.character() %>% as.numeric(),
                     predict_test_df %>% pull(response)) %>%
      performance(measure="auc") %>% .@y.values)
}

svm_start <- Sys.time()
# Store ensemble RF results as list of lists
svm_end <- lapply(training_rows, full_svm_analysis, data=eid2_model)
svm_end <- Sys.time()