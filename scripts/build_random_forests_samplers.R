# Build random forests as in build_random_forests.R using alternative resampling methods
# Note that this analysis requires significant computing power, taking approximately 60h on two Intel Xeon(R) quad-core CPUs with 36GB memory

#######################
# Load packages, data #
#######################

library(caret)
library(e1071)
library(matrixStats)
library(magrittr)
library(pROC)
library(randomForest)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)
library(stringr)
library(tibble)
library(DMwR)

# Load in data
load("data\\cov_ML_dfs_01_10_20.RData")

######################################################
# Define function to process predictive performances #
######################################################

return_sampler_values <- function(rf, validate, output = "overall"){
  
  # Predictions, excluding epidemic human spillover sequences
  predict_class_test <- Map(function(model, newdata)
    
    if (nrow(newdata)>0){
      predict(model, newdata=newdata, type="raw")
    },
    model = rf,
    newdata = validate
    
  ) %>% unlist %>% as.factor
  
  predict_prob_test <- Map(function(model, newdata)
    
    if (nrow(newdata)>0){
      predict(model, newdata = newdata, type = "prob")
    },
    model = rf,
    newdata = validate
  ) %>% bind_rows
  
  matrix_test <- confusionMatrix(predict_class_test, validate %>% bind_rows %>% pull(outcome) %>% droplevels)
  
  matrix_one_vs_all <- vector("list", length(levels(model_df$outcome)))
  for (i in seq_along(matrix_one_vs_all)) {
    positive.class <- levels(model_df$outcome)[i]
    # in the i-th iteration, use the i-th class as the positive class
    matrix_one_vs_all[[i]] <- confusionMatrix(predict_class_test, validate %>% bind_rows %>% pull(outcome) %>% droplevels, 
                                              positive = positive.class)
  }
  
  t1_df <- matrix_test$overall %>%
    round(., 3) %>%
    t() %>%
    cbind(., 
          AUC = multiclass.roc(response = validate %>% 
                                 bind_rows() %>% 
                                 pull(outcome), predictor = predict_prob_test) %>% 
            .$auc %>% 
            as.numeric() %>% 
            round(3),
          micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
          macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3)) 
  
  if (output == "classwise"){
    return(matrix_test$byClass %>% reshape2::melt())
  } else {
    return(as.data.frame(t1_df))
  }
}

##########################
# Spike protein analysis #
##########################

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model #
##################################################################################

data <- cov_spikes_df

# Set options for specific analysis
set.seed(1315)
host_cats <- c("aves", "camel", "carnivore", "human", "rodent", "swine", "yangbat", "yinbat")

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
if (!(outcome_name %in% names(data))) {
  model_df_predownsample <- data %>%
    left_join(allcov_df %>% select(childtaxa_id, !!sym(outcome_name)),
      by = c("taxid" = "childtaxa_id")
    ) %>%
    mutate(outcome = factor(!!sym(outcome_name))) %>%
    filter(!is.na(outcome))
} else {
  model_df_predownsample <- data %>%
    mutate(outcome = factor(!!sym(outcome_name))) %>%
    filter(!is.na(outcome))
}

# Exclude zoonotic and epizootic viruses
# (MERS-CoV: Middle East respiratory syndrome-related coronavirus, Betacoronavirus England 1, Human betacoronavirus 2c EMC/2012, Human betacoronavirus 2c England-Qatar/2012, Human betacoronavirus 2c Jordan-N3/2012)
# (SARS-CoV: SARS coronavirus HKU-39849, SARS coronavirus P2, Severe acute respiratory syndrome-related coronavirus)
# (SARS-CoV-2: Severe acute respiratory syndrome coronavirus 2)
# (SADS-CoV: Swine acute diarrhea syndrome coronavirus, Swine enteric alphacoronavirus, Porcine enteric alphacoronavirus GDS04, Porcine enteric alphacoronavirus)
# (other: Human enteric coronavirus strain 4408, Human enteric coronavirus 4408)

if (exclude_zoonotic == TRUE) {
  model_df %<>% filter(!(taxid %in% c(1335626, 694009, 627442, 228404, 2697049, 1263720, 1235996, 1298362, 1306931, 627439, 166124) & outcome == "human"))
  model_df %<>% filter(!(taxid %in% c(2032731, 2045491, 2018513, 2026237) & outcome == "swine"))
}

# Include nucleotide, dincucleotide and codon usage composition bias measurements, with option to exclude stop codon RSCU
if (use_stop_codons == TRUE) {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion)
} else {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias))
}

model_df %<>% filter(outcome %in% host_cats) %>% mutate(outcome = factor(outcome))

# Data thinning; for each outcome:coronavirus combination, if n > 20 downsample to 20
templist <- list() # Setup empty list

model_df_predownsample %<>% mutate(sampler = factor(paste(outcome, childtaxa_name)))

for (i in 1:nlevels(model_df_predownsample$sampler)) {
  temp <- model_df_predownsample %>%
    filter(sampler == levels(model_df_predownsample$sampler)[i])

  if (nrow(temp) > 20) {

    # Randomly sample down to 20
    temp %<>% sample_n(20)
  }

  templist[[i]] <- temp
}

model_df_datathinned <- bind_rows(templist) %>% select(-sampler)


# Create outer folds for hold-one-out validation
outer_fold_dfs_datathinned <- lapply(unique(model_df_datathinned$taxid), function(x) {
  model_df_datathinned %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices_datathinned <- lapply(outer_fold_dfs_datathinned, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})


# Retain full data for alternative resampling methods within training function
model_df <- model_df_predownsample

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})

# Specify variables used
preds <- model_df %>%
  select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>%
  remove_constant() %>%
  names()

######################
# Run random forests #
######################

timer_start <- Sys.time()

# Downsample, upsample, SMOTE
for (sampler in c("down", "up", "smote")) {

  # Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger

  # Store result as list of n models where n = number of CoV species or unranked subspecies, holding one out each time
  rf_list <- Map(function(outer, inner) {
    train(
      x = outer %>% select(preds),
      y = outer %>% pull(outcome),
      method = "ranger",
      preProc = c("center", "scale"),
      metric = "Accuracy",
      num.trees = 1000,
      importance = "impurity",
      trControl = trainControl(
        method = "repeatedcv",
        sampling = sampler,
        index = inner,
        number = 10,
        repeats = 1,
        # verboseIter = TRUE,
        classProbs = TRUE
      ),
      tuneGrid = expand.grid(
        .splitrule = "gini",

        .min.node.size = seq(from = 5, to = 20, length = 3),
        .mtry = seq(from = 5, to = 20, length = 3)
      )
    )
  },
  outer = outer_fold_dfs,
  inner = inner_fold_indices
  )

  assign(paste0("rf_list_", sampler), rf_list, envir = .GlobalEnv)
}

# No resampling
rf_list_none <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    trControl = trainControl(
      method = "repeatedcv",
      sampling = NULL,
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",

      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs,
inner = inner_fold_indices
)

# Class weighting
rf_list_classweights <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    weights = unname(case_when(
      outer$outcome == "aves" ~ (1 / 8) / table(outer$outcome)[1],
      outer$outcome == "camel" ~ (1 / 8) / table(outer$outcome)[2],
      outer$outcome == "carnivore" ~ (1 / 8) / table(outer$outcome)[3],
      outer$outcome == "human" ~ (1 / 8) / table(outer$outcome)[4],
      outer$outcome == "rodent" ~ (1 / 8) / table(outer$outcome)[5],
      outer$outcome == "swine" ~ (1 / 8) / table(outer$outcome)[6],
      outer$outcome == "yangbat" ~ (1 / 8) / table(outer$outcome)[7],
      outer$outcome == "yinbat" ~ (1 / 8) / table(outer$outcome)[8]
    )),
    trControl = trainControl(
      method = "repeatedcv",
      sampling = NULL,
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",
      # .min.node.size = 5,
      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs,
inner = inner_fold_indices
)

# Data thinning (as in build_random_forests.R)
rf_list_datathinned <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    trControl = trainControl(
      method = "repeatedcv",
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",

      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs_datathinned,
inner = inner_fold_indices_datathinned
)

#####################
# Store all outputs #
#####################

save(rf_list_datathinned,
  rf_list_none,
  rf_list_classweights,
  rf_list_up,
  rf_list_down,
  rf_list_smote,
  file = paste0("listresults_ml_spike_samplers_", format(Sys.time(), "%d_%m_%y"), ".RData")
)

#####################################
# Calculate predictive performances #
#####################################

# Specify held out data as test sets (preserving order of taxonomic ids within data-thinned and non-data-thinned datasets respectively)
validate_ordall <- lapply(unique(model_df$taxid), function(x)
  model_df_datathinned %>% filter(taxid == x)
)

validate_ordorig <- lapply(unique(model_df_datathinned$taxid), function(x)
  model_df_datathinned %>% filter(taxid == x)
)

# Produce Supplementary Table S3
bind_rows(rf_list_datathinned %>% return_sampler_values(validate = validate_ordorig),
          rf_list_none %>% return_sampler_values(validate = validate_ordall),
          rf_list_classweights %>% return_sampler_values(validate = validate_ordall),
          rf_list_down %>% return_sampler_values(validate = validate_ordall),
          rf_list_up %>% return_sampler_values(validate = validate_ordall),
          rf_list_smote_reduced %>% return_sampler_values(validate = validate_reduced)) %>% 
  cbind(set = c("data thinning", "none", "class weights", "downsample", "upsample", "SMOTE"), .) %>%
  write.csv("figures_tables//Supp_Table_3.csv")

##################################
# Whole genome sequence analysis #
##################################

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model #
##################################################################################

data <- cov_wg_df

# Set options for specific analysis
set.seed(1547)
host_cats <- c("aves", "camel", "carnivore", "human", "rodent", "swine", "yangbat", "yinbat")

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
if (!(outcome_name %in% names(data))) {
  model_df_predownsample <- data %>%
    left_join(allcov_df %>% select(childtaxa_id, !!sym(outcome_name)),
              by = c("taxid" = "childtaxa_id")
    ) %>%
    mutate(outcome = factor(!!sym(outcome_name))) %>%
    filter(!is.na(outcome))
} else {
  model_df_predownsample <- data %>%
    mutate(outcome = factor(!!sym(outcome_name))) %>%
    filter(!is.na(outcome))
}

# Exclude zoonotic and epizootic viruses
# (MERS-CoV: Middle East respiratory syndrome-related coronavirus, Betacoronavirus England 1, Human betacoronavirus 2c EMC/2012, Human betacoronavirus 2c England-Qatar/2012, Human betacoronavirus 2c Jordan-N3/2012)
# (SARS-CoV: SARS coronavirus HKU-39849, SARS coronavirus P2, Severe acute respiratory syndrome-related coronavirus)
# (SARS-CoV-2: Severe acute respiratory syndrome coronavirus 2)
# (SADS-CoV: Swine acute diarrhea syndrome coronavirus, Swine enteric alphacoronavirus, Porcine enteric alphacoronavirus GDS04, Porcine enteric alphacoronavirus)
# (other: Human enteric coronavirus strain 4408, Human enteric coronavirus 4408)

if (exclude_zoonotic == TRUE) {
  model_df %<>% filter(!(taxid %in% c(1335626, 694009, 627442, 228404, 2697049, 1263720, 1235996, 1298362, 1306931, 627439, 166124) & outcome == "human"))
  model_df %<>% filter(!(taxid %in% c(2032731, 2045491, 2018513, 2026237) & outcome == "swine"))
}

# Include nucleotide, dincucleotide and codon usage composition bias measurements, with option to exclude stop codon RSCU
if (use_stop_codons == TRUE) {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion)
} else {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, genus, taxid, childtaxa_name, accessionversion) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias))
}

model_df %<>% filter(outcome %in% host_cats) %>% mutate(outcome = factor(outcome))

# Data thinning; for each outcome:coronavirus combination, if n > 20 downsample to 20
templist <- list() # Setup empty list

model_df_predownsample %<>% mutate(sampler = factor(paste(outcome, childtaxa_name)))

for (i in 1:nlevels(model_df_predownsample$sampler)) {
  temp <- model_df_predownsample %>%
    filter(sampler == levels(model_df_predownsample$sampler)[i])
  
  if (nrow(temp) > 20) {
    
    # Randomly sample down to 20
    temp %<>% sample_n(20)
  }
  
  templist[[i]] <- temp
}

model_df_datathinned <- bind_rows(templist) %>% select(-sampler)


# Create outer folds for hold-one-out validation
outer_fold_dfs_datathinned <- lapply(unique(model_df_datathinned$taxid), function(x) {
  model_df_datathinned %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices_datathinned <- lapply(outer_fold_dfs_datathinned, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})


# Retain full data for alternative resampling methods within training function
model_df <- model_df_predownsample

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})

# Specify variables used
preds <- model_df %>%
  select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>%
  remove_constant() %>%
  names()

######################
# Run random forests #
######################

timer_start <- Sys.time()

# Downsample, upsample, SMOTE
for (sampler in c("down", "up", "smote")) {

  # Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger

  # Store result as list of n models where n = number of CoV species or unranked subspecies, holding one out each time
  rf_list <- Map(function(outer, inner) {
    train(
      x = outer %>% select(preds),
      y = outer %>% pull(outcome),
      method = "ranger",
      preProc = c("center", "scale"),
      metric = "Accuracy",
      num.trees = 1000,
      importance = "impurity",
      trControl = trainControl(
        method = "repeatedcv",
        sampling = sampler,
        index = inner,
        number = 10,
        repeats = 1,
        # verboseIter = TRUE,
        classProbs = TRUE
      ),
      tuneGrid = expand.grid(
        .splitrule = "gini",

        .min.node.size = seq(from = 5, to = 20, length = 3),
        .mtry = seq(from = 5, to = 20, length = 3)
      )
    )
  },
  outer = outer_fold_dfs,
  inner = inner_fold_indices
  )

  assign(paste0("rf_list_", sampler), rf_list, envir = .GlobalEnv)
}

# No resampling
rf_list_none <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    trControl = trainControl(
      method = "repeatedcv",
      sampling = NULL,
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",

      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs,
inner = inner_fold_indices
)

# Class weighting
rf_list_classweights <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    weights = unname(case_when(
      outer$outcome == "aves" ~ (1 / 8) / table(outer$outcome)[1],
      outer$outcome == "camel" ~ (1 / 8) / table(outer$outcome)[2],
      outer$outcome == "carnivore" ~ (1 / 8) / table(outer$outcome)[3],
      outer$outcome == "human" ~ (1 / 8) / table(outer$outcome)[4],
      outer$outcome == "rodent" ~ (1 / 8) / table(outer$outcome)[5],
      outer$outcome == "swine" ~ (1 / 8) / table(outer$outcome)[6],
      outer$outcome == "yangbat" ~ (1 / 8) / table(outer$outcome)[7],
      outer$outcome == "yinbat" ~ (1 / 8) / table(outer$outcome)[8]
    )),
    trControl = trainControl(
      method = "repeatedcv",
      sampling = NULL,
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",
      # .min.node.size = 5,
      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs,
inner = inner_fold_indices
)

# Data thinning (as in build_random_forests.R)
rf_list_datathinned <- Map(function(outer, inner) {
  train(
    x = outer %>% select(preds),
    y = outer %>% pull(outcome),
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Accuracy",
    num.trees = 1000,
    importance = "impurity",
    trControl = trainControl(
      method = "repeatedcv",
      index = inner,
      number = 10,
      repeats = 1,
      # verboseIter = TRUE,
      classProbs = TRUE
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",

      .min.node.size = seq(from = 5, to = 20, length = 3),
      .mtry = seq(from = 5, to = 20, length = 3)
    )
  )
},
outer = outer_fold_dfs_datathinned,
inner = inner_fold_indices_datathinned
)

#####################
# Store all outputs #
#####################

save(rf_list_datathinned,
  rf_list_none,
  rf_list_classweights,
  rf_list_up,
  rf_list_down,
  rf_list_smote,
  file = paste0("listresults_ml_wg_samplers_", format(Sys.time(), "%d_%m_%y"), ".RData")
)

#####################################
# Calculate predictive performances #
#####################################

# Specify held out data as test sets (preserving order of taxonomic ids within data-thinned and non-data-thinned datasets respectively)
validate_ordall <- lapply(unique(model_df$taxid), function(x)
  model_df_datathinned %>% filter(taxid == x)
)

validate_ordorig <- lapply(unique(model_df_datathinned$taxid), function(x)
  model_df_datathinned %>% filter(taxid == x)
)

# Produce Supplementary Table S4
bind_rows(rf_list_datathinned %>% return_sampler_values(validate = validate_ordorig),
          rf_list_none %>% return_sampler_values(validate = validate_ordall),
          rf_list_classweights %>% return_sampler_values(validate = validate_ordall),
          rf_list_down %>% return_sampler_values(validate = validate_ordall),
          rf_list_up %>% return_sampler_values(validate = validate_ordall),
          rf_list_smote_reduced %>% return_sampler_values(validate = validate_reduced)) %>% 
  cbind(set = c("data thinning", "none", "class weights", "downsample", "upsample", "SMOTE"), .) %>%
  write.csv("figures_tables//Supp_Table_4.csv")