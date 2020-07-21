library("DALEX")

# pull out best model
rf_list[[1]] %>% .$tuned_obj %>% .$best.model

# define predict function to extract human prediction
pred_rf_human <- function(model, data){
  predict(model, newdata=data, type="prob") %>%
    as.data.frame %>%
    pull(human)
}

# Calculate variable attribution
bd_rf <- variable_attribution(explain(model = rf_list[[1]] %>% .$tuned_obj %>% .$best.model,
                                      data = rf_list[[1]] %>% .$train %>% select(-outcome),
                                      y = rf_list[[1]] %>% .$train %>% pull(outcome),
                                      predict_function = pred_rf_human,
                                      label = "random forest"),
                              new_observation = rf_list[[1]] %>% .$train %>% slice(1) %>% select(-outcome),
                              type = "break_down_interactions")

plot(bd_rf, max_features = 10)

# # Include strands for all data points
# bd_rf_distr <- variable_attribution(explain_rf_v6,
#                                     new_observation = henry, 
#                                     type = "break_down",
#                                     keep_distributions = TRUE)
# 
# plot(bd_rf_distr, plot_distributions = TRUE) 


#### PRACTICE WITH VIRULENCE PAPER RF MODEL

setwd(r"(C:\Users\Liam\PHD\BACKUPS\2 Oct 2015 Full Working\PHD\Writeups\Virulence chapter\Ch1_virulence_backup_27_8_19\Ch1 Virulence analysis resubmission)")

# Read in and process data, selecting only those viruses with reliable virulence data for includsion in analysis
vir <- read.csv("Virulence_data_Brierley_et_al_21_6_19.csv", sep=",", header=T, stringsAsFactors = TRUE)
vir <- subset(vir, Included==1)

# Remove 'Unassigned' as a family and rearrange by genome type
vir$Family <- droplevels(vir$Family)
vir$Family <- factor(vir$Family, levels=c("Arenaviridae","Bornaviridae","Filoviridae","Hantaviridae","Nairoviridae","Orthomyxoviridae","Paramyxoviridae","Peribunyaviridae","Phenuiviridae","Pneumoviridae","Rhabdoviridae","Astroviridae","Caliciviridae","Coronaviridae","Flaviviridae","Hepeviridae","Picornaviridae","Togaviridae","Picobirnaviridae","Reoviridae","Retroviridae"), labels=c("Arena-","Borna-","Filo-","Hanta-","Nairo-","Orthomyxo-","Paramyxo-","Peribunya-","Phenui-","Pneumo-","Rhabdo-","Astro-","Calici-","Corona-","Flavi-","Hepe-","Picorna-","Toga-","Picobirna-","Reo-","Retro-"))

# Convert 0/1 columns to factors
factor_columns <- c(names(vir[,grepl( "^Tr\\.|^Fatal\\.|^Tp\\.|^Host\\." , names(vir))]))
vir[factor_columns] <- lapply(vir[factor_columns], factor)

# Create True Skill Statistic function
trueskillstat <- function(matrix){
  as.numeric(matrix$byClass[1] + matrix$byClass[2] - 1)
}

# Specify predictor variables to use in models as a formula
formula_used <- formula(Virulence ~ Family+Genome+
                          Tr.level+
                          Tr.primary+
                          Tr.direct.contact+Tr.faecal.oral+Tr.respiratory+Tr.vector+Tr.multiple+Tr.food.borne+Tr.vertical+
                          Tp.primary+
                          Tp.vascular+Tp.circulatory+Tp.gastrointestinal+Tp.hepatic+Tp.neural+Tp.respiratory+Tp.cardiac+Tp.joints+Tp.renal+Tp.reproductive+Tp.sensory+Tp.skin+Tp.muscular+Tp.endocrine+Tp.multiple+
                          Host.range+Host.human.only+Host.nh.primates+Host.other.mammal+Host.bird)

# Set random seed for reproducible training/test partition
set.seed(1353)

# Establish training and test partition based on severity and taxonomy at 3:1 ratio (if < 4 viruses exist for a given virulence-genus combination, they default to the training set)
training_row_nos = suppressWarnings(createDataPartition(paste(vir$Virulence, vir$Genus), p = 0.75, list = FALSE))
vir_train <- vir[training_row_nos,]
vir_test <- vir[-training_row_nos,]
# Store test set virus species names
vir_test_names <- vir_test$Species

# Impute missing data, itself from a random forest
invisible(capture.output(data_imp_train <- missForest(vir_train[, all.vars(formula_used)])))
vir_train <- data_imp_train$ximp
invisible(capture.output(data_imp_test <- missForest(vir_test[, all.vars(formula_used)])))
vir_test <- data_imp_test$ximp

# Build random forest model from 5000 individual bootstrapped trees
random_forest <- randomForest(formula_used, 
                              dat= vir_train,
                              importance = TRUE, 
                              na.action = na.exclude, 
                              keep.forest = TRUE, 
                              ntree = 5000, 
                              proximity = TRUE,
                              nodesize = 5)

# define predict function to extract severe prediction
pred_rf_severe <- function(model, data){
  predict(model, newdata=data, type="prob") %>%
    as.data.frame %>%
    pull(severe)
}

bd_usutu <- variable_attribution(explain(model = random_forest,
                                      data = vir_train %>% select(labels(terms(formula_used))),
                                      y = vir_train %>% pull(Virulence),
                                      predict_function = pred_rf_severe,
                                      label = "Brierley et al. random forest, profile: Usutu virus"),
                              new_observation = vir %>% filter(Species == "Usutu virus") %>% select(labels(terms(formula_used))),
                              type = "break_down")

plot(bd_usutu, max_features = 10)
