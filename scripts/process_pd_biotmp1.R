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

load("listresults_ml_vector_spike_17_07_20.RData")

vars_to_pd <- lapply(rf_list, function(x)
  varImp(x)$importance %>% rownames_to_column("name") %>% mutate(Overall = Overall/100) %>% rename(relGini = Overall)
) %>% 
  bind_rows() %>%
  group_by(name) %>% 
  summarise(mean = mean(relGini)) %>% 
  arrange(-mean) %>%
  slice(1:5) %>%
  pull(name)

PD_extractor <- function(model){
  
  temp_list <- vector("list", length(vars_to_pd)) %>% setNames(vars_to_pd)
  
  for (i in 1:length(vars_to_pd)){
    
    list_of_classes <- lapply(levels(model_df$outcome), function(y) 
      pdp::partial(model, 
                   pred.var = vars_to_pd[i], 
                   pred.grid = seq(from = min(model_df %>% pull(vars_to_pd[i])), 
                                   to = max(model_df %>% pull(vars_to_pd[i])), 
                                   length.out = 30) %>% 
                     as.data.frame() %>% 
                     setNames(vars_to_pd[i]),
                   which.class = y,
                   prob = TRUE,
                   progress = progress_text(char = "-"))) %>%
      setNames(levels(rf_list$Resample1$train$outcome))
    
    pred_by_class <- join_all(list_of_classes, by=vars_to_pd[i], type='left') %>% set_colnames(c(vars_to_pd[i], levels(model_df$outcome))) 
    
    temp_list[[i]] <- pred_by_class
  }
  
  return(temp_list)
}

timer_start <- Sys.time()
list_PD <- lapply(rf_list[1:2], function(x) PD_extractor(x))
timer_end <- Sys.time()

save(timer_start, timer_end, list_PD,
     file=paste0("listresults_pd_spike_", format(Sys.time(), "%d_%m_%y"), ".RData"))