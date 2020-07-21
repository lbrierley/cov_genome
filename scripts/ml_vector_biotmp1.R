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
library(ranger)

# Load in previous ML results
load("cov_ML_dfs_28_05_20.RData")

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model # Spike protein
##################################################################################

data <- cov_spikes_df

# Set options
set.seed(1315)
outcome_name <- "group_name"
use_stop_codons <- TRUE
min_n_seq_category <- 60

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
if (!(outcome_name %in% names(data))){
  model_df_predownsample <- data %>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome_name)),
                                               by = c("taxid" = "childtaxa_id")) %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
} else {
  model_df_predownsample <- data %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
}

# Downsample outcome:CoV species
# For each outcome:CoV species combination, if n > 20 downsample to 20
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

# Specify variables used
preds <- model_df %>% select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>% remove_constant %>% names

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x)
  model_df %>% filter(taxid != x)
)

holdout_dfs <- lapply(unique(model_df$taxid), function(x)
  model_df %>% filter(taxid == x)
)

holdout_dfs_predownsample <- lapply(unique(model_df$taxid), function(x)
  model_df_predownsample %>% filter(outcome %in% unique(model_df$outcome) & taxid == x)
)

set.seed(141)
holdout_dfs_minimal <- lapply(unique(model_df$taxid), function(x)
  model_df_predownsample %>% filter(outcome %in% unique(model_df$outcome) & taxid == x) %>% group_by(outcome) %>% sample_n(1)
)

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x)
  createMultiFolds(x$outcome, k = 10, times = 1)
)

######################
# Run random forests #
######################

# Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger

timer_start <- Sys.time()

# Store result as list of n ensemble models where n = number of CoV species, holding one out each time
rf_list <- Map(function(outer, inner) 
  
  train(x = outer %>% select(preds),
        y = outer %>% pull(outcome),
        method = "ranger",
        preProc = c("center", "scale"),
        metric = "Accuracy",
        num.trees = 1000,
        importance = "impurity",
        trControl = trainControl(method = "repeatedcv", 
                                 index = inner,
                                 number = 10,
                                 repeats = 1,
                                 #verboseIter = TRUE,
                                 classProbs = TRUE),
        tuneGrid = expand.grid(
          .splitrule = "gini",
          #.min.node.size = 5,
          .min.node.size = seq(from = 5, to = 20, length = 3),
          .mtry = seq(from = 5, to = 20, length = 3))
  ),
  outer = outer_fold_dfs,
  inner = inner_fold_indices
)

timer_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df, model_df_predownsample, use_stop_codons,
     #svm_start, svm_end, svm_list,
     timer_start, timer_end, rf_list,
     file=paste0("listresults_ml_vector_spike_", format(Sys.time(), "%d_%m_%y"), ".RData"))




##################################################################################
# Prepare data frame for modelling and define set of variables used in the model # Wgs
##################################################################################

data <- cov_wg_df

# Set options
set.seed(1547)
outcome_name <- "group_name"
hold_out_species <- TRUE
use_stop_codons <- TRUE
min_n_seq_category <- 40

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
if (!(outcome_name %in% names(data))){
  model_df_predownsample <- data %>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome_name)),
                                               by = c("taxid" = "childtaxa_id")) %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
} else {
  model_df_predownsample <- data %>%
    mutate(outcome = factor(!! sym(outcome_name))) %>%
    filter(!is.na(outcome))
}

# Downsample outcome:CoV species
# For each outcome:CoV species combination, if n > 20 downsample to 20
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

# Specify variables used
preds <- model_df %>% select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>% remove_constant %>% names

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x)
  model_df %>% filter(taxid != x)
)

holdout_dfs <- lapply(unique(model_df$taxid), function(x)
  model_df %>% filter(taxid == x)
)

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x)
  createMultiFolds(x$outcome, k = 10, times = 1)
)

######################
# Run random forests #
######################

# Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger

timer_start <- Sys.time()

# Store result as list of n ensemble models where n = number of CoV species, holding one out each time
rf_list <- Map(function(outer, inner) 
  
  train(x = outer %>% select(preds),
        y = outer %>% pull(outcome),
        method = "ranger",
        preProc = c("center", "scale"),
        metric = "Accuracy",
        num.trees = 1000,
        importance = "impurity",
        trControl = trainControl(method = "repeatedcv", 
                                 index = inner,
                                 number = 10,
                                 repeats = 1,
                                 #verboseIter = TRUE,
                                 classProbs = TRUE),
        tuneGrid = expand.grid(
          .splitrule = "gini",
          #.min.node.size = 5,
          .min.node.size = seq(from = 5, to = 20, length = 3),
          .mtry = seq(from = 5, to = 20, length = 3))
  ),
  outer = outer_fold_dfs,
  inner = inner_fold_indices
)

timer_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df, model_df_predownsample, use_stop_codons,
     #svm_start, svm_end, svm_list,
     timer_start, timer_end, rf_list,
     file=paste0("listresults_ml_vector_wg_", format(Sys.time(), "%d_%m_%y"), ".RData"))