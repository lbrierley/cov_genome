# Note that this analysis requires significant computing power, taking approximately 10h on two Intel Xeon(R) quad-core CPUs with 36GB memory
# Saved analytical outputs are provided in the 'outputs' folder

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

# Load in data
load("data\\cov_ML_dfs_01_10_20.RData")

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

model_df <- bind_rows(templist) %>% select(-sampler)

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

# Specify variables used
preds <- model_df %>%
  select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>%
  remove_constant() %>%
  names()

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})

######################
# Run random forests #
######################

# Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger
timer_start <- Sys.time()

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

timer_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df, model_df_predownsample, use_stop_codons,
  timer_start, timer_end, rf_list,
  file = paste0("outputs\\listresults_ml_spike_", format(Sys.time(), "%d_%m_%y"), ".RData")
)

# Calculate and store variable importance
varimp_order_spike <- lapply(rf_list, function(x) {
  varImp(x)$importance %>%
    rownames_to_column("name") %>%
    mutate(Overall = Overall / 100) %>%
    rename(relGini = Overall)
}) %>%
  bind_rows() %>%
  group_by(name) %>%
  summarise(mean = mean(relGini)) %>%
  arrange(-mean) %>%
  pull(name)

saveRDS(varimp_order_spike, "outputs\\varimp_order_spike.rds")


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

model_df <- bind_rows(templist) %>% select(-sampler)

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

# Specify variables used
preds <- model_df %>%
  select(-outcome, -taxid, -childtaxa_name, -accessionversion, -genus) %>%
  remove_constant() %>%
  names()

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid != x)
})

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x) {
  createMultiFolds(x$outcome, k = 10, times = 1)
})

######################
# Run random forests #
######################

# Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger
timer_start <- Sys.time()

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

timer_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df, model_df_predownsample, use_stop_codons,
  timer_start, timer_end, rf_list,
  file = paste0("outputs\\listresults_ml_wg_", format(Sys.time(), "%d_%m_%y"), ".RData")
)

# Calculate and store variable importance
varimp_order_wg <- lapply(rf_list, function(x) {
  varImp(x)$importance %>%
    rownames_to_column("name") %>%
    mutate(Overall = Overall / 100) %>%
    rename(relGini = Overall)
}) %>%
  bind_rows() %>%
  group_by(name) %>%
  summarise(mean = mean(relGini)) %>%
  arrange(-mean) %>%
  pull(name)

saveRDS(varimp_order_wg, "outputs\\varimp_order_wg.rds")
