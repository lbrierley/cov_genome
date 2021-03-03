# Note that this analysis requires significant computing power, taking approximately 20-30h on two Intel Xeon(R) quad-core CPUs with 36GB memory
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
library(ranger)

# Define partial dependence extractor function

PD_extractor <- function(model) {
  temp_list <- vector("list", length(vars_to_pd)) %>% setNames(vars_to_pd)

  for (i in 1:length(vars_to_pd)) {
    list_of_classes <- lapply(levels(model_df$outcome), function(y) {
      pdp::partial(model,
        pred.var = vars_to_pd[i],
        pred.grid = seq(
          from = min(model_df %>% pull(vars_to_pd[i])),
          to = max(model_df %>% pull(vars_to_pd[i])),
          length.out = 30
        ) %>%
          as.data.frame() %>%
          setNames(vars_to_pd[i]),
        which.class = y,
        prob = TRUE
      )
    })

    pred_by_class <- plyr::join_all(list_of_classes, by = vars_to_pd[i], type = "left") %>% set_colnames(c(vars_to_pd[i], levels(model_df$outcome)))

    temp_list[[i]] <- pred_by_class
  }

  return(temp_list)
}

# Load and set list of features to calculate partial dependence for (four most important features by mean Gini impurity)

varimp_order_spike <- readRDS("outputs\\varimp_order_spike.rds")
varimp_order_wg <- readRDS("outputs\\varimp_order_wg.rds")

vars_to_pd <- unique(c(varimp_order_spike[1:4], varimp_order_wg[1:4]))

# Calculate and store partial dependence profiles, spike proteins

load("outputs\\listresults_ml_spike_02_10_20.RData")

timer_start <- Sys.time()
list_PD <- lapply(rf_list, function(x) PD_extractor(x))
timer_end <- Sys.time()

save(timer_start, timer_end, list_PD,
  file = paste0("outputs\\listresults_pd_spike_", format(Sys.time(), "%d_%m_%y"), ".RData")
)

# Calculate and store partial dependence profiles,  whole genomes

load("outputs\\outputslistresults_ml_wg_02_10_20.RData")

timer_start <- Sys.time()
list_PD <- lapply(rf_list, function(x) PD_extractor(x))
timer_end <- Sys.time()

save(timer_start, timer_end, list_PD,
  file = paste0("outputs\\listresults_pd_wg_", format(Sys.time(), "%d_%m_%y"), ".RData")
)
