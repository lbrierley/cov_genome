#################
# Load packages #
#################

library(caret)
library(e1071)
library(ggalluvial)
library(matrixStats)
library(patchwork)
library(pROC)
library(randomForest)
library(aricode)

#######################################
# Load elements common to all results #
#######################################

# Read in list of most important variables for use in PD
vars_to_pd_spike <- readRDS("outputs\\varimp_order_spike.rds")[1:4]
vars_to_pd_wg <- readRDS("outputs\\varimp_order_wg.rds")[1:4]

# Read in tax matcher to access metadata raw host strings
load(file = "data\\tax_matcher.RData")

#######################################
# Load spike protein analysis results #
#######################################

load("outputs\\listresults_ml_spike_02_10_20.RData")

# Set unclassified genus
model_df %<>% replace_na(list(genus = "unclassified"))

# Extract parameter optimisation
gridsearch <- lapply(rf_list, function(x) {
  x$results
}) %>%
  bind_rows() %>%
  mutate(
    min.node.size = factor(min.node.size, levels = c("5", "12.5", "20")),
    mtry = factor(mtry, levels = c("5", "12.5", "20"))
  )

# Select holdout used
validate_used <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid == x)
})

# Extract class predictions (excluding zoonotic sequences)
predict_class_test <- Map(function(model, newdata) {
  if (nrow(newdata) > 0) {
    levels(model_df$outcome)[predict(model, newdata = newdata, type = "raw")]
  }
},
model = rf_list,
newdata = validate_used
) %>%
  unlist() %>%
  as.factor()

# Extract probabilistic predictions (excluding zoonotic sequences)
predict_prob_test <- Map(function(model, newdata) {
  if (nrow(newdata) > 0) {
    predict(model, newdata = newdata, type = "prob")
  }
},
model = rf_list,
newdata = validate_used
) %>% bind_rows()

# Calculate confusion matrix
matrix_test <- confusionMatrix(predict_class_test, validate_used %>% bind_rows() %>% pull(outcome) %>% droplevels())

# Calculate one-vs-all confusion matrices
# Sourced from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/ (Matthias Döring, 2018)
matrix_one_vs_all <- vector("list", length(levels(model_df$outcome)))
for (i in seq_along(matrix_one_vs_all)) {
  positive.class <- levels(model_df$outcome)[i]
  # in the i-th iteration, use the i-th class as the positive class
  matrix_one_vs_all[[i]] <- confusionMatrix(predict_class_test, validate_used %>% bind_rows() %>% pull(outcome) %>% droplevels(),
    positive = positive.class
  )
}


# Collect equivalent test set for zoonotic sequences (i.e. MERS-CoV, SARS-CoV, SARS-CoV-2) sampled from humans
set.seed(1223)
validate_zoonotic <- model_df_predownsample %>%
  filter(taxid %in% c(1335626, 694009, 627442, 228404, 2697049, 1263720, 1235996, 1298362, 1306931, 627439, 166124) & outcome == "human") %>%
  split(.$taxid) %>%
  lapply(function(x) {
    x %>% sample_n(ifelse(nrow(x) < 20, nrow(x), 20))
  }) %>%
  bind_rows() %>%
  mutate(outcome = factor(outcome, levels = levels(model_df$outcome)))

# Extract probabilistic predictions for zoonotic sequences (grand mean probabilities over all random forests)
predict_prob_zoonotic <- lapply(rf_list, function(model) {
  validate_zoonotic %>%
    select(taxid, accessionversion) %>%
    cbind(predict(model, newdata = validate_zoonotic, type = "prob"))
}) %>%
  bind_rows() %>%
  group_by(taxid, accessionversion) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(flag = case_when(
    taxid %in% c(1335626, 1263720, 1235996, 1298362, 1306931) ~ "MERS-CoV",
    taxid == 2697049 ~ "SARS-CoV-2",
    taxid %in% c(694009, 627442, 228404) ~ "SARS-CoV"
  ))

# Calculate variable importance
varimp_spike <- lapply(rf_list, function(x) {
  varImp(x)$importance %>%
    rownames_to_column("name") %>%
    mutate(Overall = Overall / 100) %>%
    rename(relGini = Overall)
}) %>%
  bind_rows() %>%
  mutate(name = str_replace_all(name, c(
    "_Bias" = "",
    "T" = "U",
    "_p" = " \\(p",
    "1$" = "1\\)",
    "2$" = "2\\)",
    "3$" = "3\\)"
  )))

varimp_order_spike <- varimp_spike %>%
  group_by(name) %>%
  summarise(mean = mean(relGini)) %>%
  mutate(name = fct_reorder(name, -mean)) %>%
  pull(name) %>%
  levels()

# Read in partial dependence profiles
load("outputs\\listresults_pd_spike_03_10_20.RData")

##################################
# Create/save tables and figures #
##################################

# Tables

# Table 1, descriptive stats - mean, SD over each coronavirus species and unranked subspecies
enc_per_spp <- model_df_predownsample %>%
  filter(accessionversion %in% model_df$accessionversion) %>%
  group_by(genus, accessionversion) %>%
  summarise(spp_mean = mean(enc)) %>%
  ungroup()

bind_rows(
  enc_per_spp %>%
    group_by(genus) %>%
    summarise(
      mean = mean(spp_mean),
      sd = sd(spp_mean)
    ),
  data.frame(genus = "Total", mean = mean(enc_per_spp$spp_mean), sd = sd(enc_per_spp$spp_mean))
) %>%
  write.csv("figures_tables\\Table_1_spike.csv")

# Table 2, predictive performance of random forest models
matrix_test$overall %>%
  round(., 3) %>%
  t() %>%
  cbind(.,
    AUC = multiclass.roc(response = validate_used %>%
      bind_rows() %>%
      pull(outcome), predictor = predict_prob_test) %>%
      .$auc %>%
      as.numeric() %>%
      round(3),
    micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
    macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3)
  ) %>%
  write.csv("figures_tables\\Table_2_spike.csv")

# Supplementary Tables

# Supplementary Table S2, number of sequences per host category
model_df %>%
  group_by(outcome) %>%
  summarise(n_sequences = n(), n_species = n_distinct(taxid)) %>%
  write.csv("figures_tables\\Supp_Table_2_spike.csv")

# Supplementary Table 6, precision/recall by class
matrix_test$byClass %>% write.csv("figures_tables\\Supp_Table_6.csv")

# Figures

# Figure 1, heatmap
png("figures_tables\\Figure_1A.png", width = 12, height = 9, units = "in", res = 320)
heatmap_spike <- model_df %>%
  as.data.frame() %>%
  set_rownames(.$accessionversion) %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][A|C|G|U]$")) %>%
  janitor::remove_constant() %>%
  as.matrix() %>%
  t() %>%
  heatmap.2(
    density.info = "none",
    trace = "none",
    margins = c(1, 5),
    dendrogram = "col",
    Rowv = "NA",
    labCol = NA,
    ColSideColors = model_df %>% mutate(colsidecol = case_when( # Set side colours using same genus colours as ggplots elsewhere
      genus == "Alphacoronavirus" ~ "#D55E00",
      genus == "Betacoronavirus" ~ "#F0E442",
      genus == "Deltacoronavirus" ~ "#0072B2",
      genus == "Gammacoronavirus" ~ "#009E73",
      genus == "unclassified" ~ "#999999",
    )) %>% pull(colsidecol),
    cexRow = 1.1,
    col = colorRampPalette(c("dodgerblue", "gray10", "firebrick2"))(n = 65),
    breaks = c(seq(0, 0.95, length = 25), seq(0.951, 1.05, length = 6), seq(1.051, 3.7, length = 35)),
    keysize = 0.75, key.par = list(cex = 0.5), key.title = NA, key.xlab = "RSCU",
    lhei = c(2, 10), lwid = c(2, 10)
  )
par(lend = 1) # square line ends for the color legend
legend(
  cex = 0.8, x = -0.05, y = 0.9, xpd = TRUE,
  legend = model_df$genus %>% unique() %>% sort(), fill = c("#D55E00", "#F0E442", "#0072B2", "#009E73", "#999999"), ncol = 1
)
dev.off()

# Compare dendogram clustering to genus using Normalised Mutal Information
heatmap_spike %>%
  .$colDendrogram %>%
  as.hclust() %>%
  cutree(k = 10) %>%
  as.data.frame() %>%
  set_colnames("cluster") %>%
  rownames_to_column("accessionversion") %>%
  left_join(model_df %>% select(accessionversion, genus)) %>%
  filter(genus != "unclassified") %>%
  with(., NMI(cluster, genus))

# Figure 2, model predictions for held-out coronavirus sequences
valid_df_raw <- data.frame(
  host_category = validate_used %>% bind_rows() %>% pull(outcome),
  accessionversion = validate_used %>% bind_rows() %>% pull(accessionversion),
  genus = validate_used %>% bind_rows() %>% pull(genus),
  childtaxa_name = validate_used %>% bind_rows() %>% pull(childtaxa_name),
  predict_prob_test
)

valid_df_ordered_spike <- rbind(
  valid_df_raw %>% filter(host_category == "aves") %>% arrange(-aves),
  valid_df_raw %>% filter(host_category == "camel") %>% arrange(-camel),
  valid_df_raw %>% filter(host_category == "carnivore") %>% arrange(-carnivore),
  valid_df_raw %>% filter(host_category == "human") %>% arrange(factor(childtaxa_name, levels = valid_df_raw %>%
    filter(host_category == "human") %>%
    arrange(-human) %>%
    pull(childtaxa_name) %>%
    unique()), -human),
  valid_df_raw %>% filter(host_category == "rodent") %>% arrange(-rodent),
  valid_df_raw %>% filter(host_category == "swine") %>% arrange(-swine),
  valid_df_raw %>% filter(host_category == "yangbat") %>% arrange(-yangbat),
  valid_df_raw %>% filter(host_category == "yinbat") %>% arrange(-yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  select(-genus) %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))

F2_spike <- valid_df_ordered_spike %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% mean(), y = 0.05, label = "HKU1", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% mean(), y = 0.05, label = "NL63", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% mean(), y = 0.05, label = "229E", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_spike %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% mean(), y = 0.05, label = "OC43", fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Figure 3, model predictions for zoonotic coronavirus sequences
F3_spike <- rbind(
  predict_prob_zoonotic %>% filter(flag == "MERS-CoV") %>% arrange(-camel),
  predict_prob_zoonotic %>% filter(flag == "SARS-CoV") %>% arrange(-rodent),
  predict_prob_zoonotic %>% filter(flag == "SARS-CoV-2") %>% arrange(-yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  melt(id.vars = c("accessionversion", "taxid", "flag")) %>%
  mutate_at(vars(variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  ))) %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_grid(~flag, scales = "free_x", space = "free_x")

# Supplementary Figures

# Supplementary Figure S2 - grid search results
SF2_spike <- ggplot(gridsearch, aes(y = Accuracy, x = min.node.size, fill = mtry)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Cross-validation (inner loop) accuracy")

# Supplementary Figure S3, dinucleotide bias values by position
SF3_spike <- model_df %>%
  rename_with(., ~ str_replace_all(., c(
    "_Bias" = "",
    "T" = "U",
    "_p" = " \\(p",
    "1$" = "1\\)",
    "2$" = "2\\)",
    "3$" = "3\\)"
  ))) %>%
  select(matches("^[A|C|G|U][A|C|G|U] ")) %>%
  melt() %>%
  mutate(
    dinuc = substr(variable, 1, 2),
    position = substr(variable, 6, 6),
    position = case_when(
      position == "1" ~ "1-2",
      position == "2" ~ "2-3",
      position == "3" ~ "3-1"
    )
  ) %>%
  ggplot(aes(y = value, x = position, fill = position)) +
  geom_boxplot(outlier.shape = NA, width = 0.9) +
  facet_wrap(. ~ variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  ylab("Dinucleotide bias") +
  geom_hline(yintercept = 1, alpha = 0.3, color = "grey50", lty = "dashed", size = 1) +
  facet_wrap(~dinuc)

# Supplementary Figure S4, model predictions for held-out coronavirus sequences blocked by genus
valid_df_ordered_spike_genus <- rbind(
  valid_df_raw %>% filter(host_category == "aves") %>% arrange(genus, -aves),
  valid_df_raw %>% filter(host_category == "camel") %>% arrange(genus, -camel),
  valid_df_raw %>% filter(host_category == "carnivore") %>% arrange(genus, -carnivore),
  valid_df_raw %>% filter(host_category == "human") %>% arrange(genus, -human),
  valid_df_raw %>% filter(host_category == "rodent") %>% arrange(genus, -rodent),
  valid_df_raw %>% filter(host_category == "swine") %>% arrange(genus, -swine),
  valid_df_raw %>% filter(host_category == "yangbat") %>% arrange(genus, -yangbat),
  valid_df_raw %>% filter(host_category == "yinbat") %>% arrange(genus, -yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index", "genus")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))


known_combs <- valid_df_ordered_spike_genus %>% distinct(host_category, genus)
spike_genus_annots <- list()
for (i in 1:nrow(known_combs)) {
  spike_genus_annots[[i]] <- data.frame(
    host_category = known_combs$host_category[i],
    genus = known_combs$genus[i],
    xmin = valid_df_ordered_spike_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% min() - 0.5,
    xmax = valid_df_ordered_spike_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% max() + 0.5,
    xmean = valid_df_ordered_spike_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% mean()
  ) %>%
    mutate(label = case_when(
      genus == "Alphacoronavirus" ~ "alpha",
      genus == "Betacoronavirus" ~ "beta",
      genus == "Gammacoronavirus" ~ "gamma",
      genus == "Deltacoronavirus" ~ "delta",
      genus == "unclassified" ~ "U"
    ))
}

spike_genus_annots <- bind_rows(spike_genus_annots)

SF4_spike_genus <- valid_df_ordered_spike_genus %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = spike_genus_annots,
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = xmin, xmax = xmax, ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = spike_genus_annots, size = 3,
    aes(x = xmean, y = 0.05, label = label, fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Supplementary Figure S5, model predictions for held-out coronavirus sequences blocked by species
valid_df_ordered_spike_spp <- rbind(
  rearrange_to_plot_pred_fig_spp("aves"),
  rearrange_to_plot_pred_fig_spp("camel"),
  rearrange_to_plot_pred_fig_spp("carnivore"),
  rearrange_to_plot_pred_fig_spp("human"),
  rearrange_to_plot_pred_fig_spp("rodent"),
  rearrange_to_plot_pred_fig_spp("swine"),
  rearrange_to_plot_pred_fig_spp("yangbat"),
  rearrange_to_plot_pred_fig_spp("yinbat")
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  select(-genus) %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))


known_combs <- valid_df_ordered_spike_spp %>% distinct(host_category, childtaxa_name)
spike_spp_annots <- list()
for (i in 1:nrow(known_combs)) {
  spike_spp_annots[[i]] <- data.frame(
    host_category = known_combs$host_category[i],
    childtaxa_name = known_combs$childtaxa_name[i],
    xmin = valid_df_ordered_spike_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% min() - 0.5,
    xmax = valid_df_ordered_spike_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% max() + 0.5,
    xmean = valid_df_ordered_spike_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% mean()
  )
}

spike_spp_annots <- bind_rows(spike_spp_annots) %>%
  mutate(childtaxa_name = factor(childtaxa_name, levels = unique(childtaxa_name))) %>%
  group_by(childtaxa_name) %>%
  mutate(label = cur_group_id()) %>%
  ungroup()

SF5_spike_spp <- valid_df_ordered_spike_spp %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = spike_spp_annots,
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = xmin, xmax = xmax, ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.3
  ) +
  geom_text(
    data = spike_spp_annots, size = 1.4, angle = 90,
    aes(x = xmean, y = -0.05, label = label, fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0.05)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Supplementary Figure S6 & S7, partial dependence
for (i in 1:length(vars_to_pd_spike)) {
  plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_spike[i]]]) %>%
    bind_rows() %>%
    rename_with(., ~ str_replace_all(., c(
      "_Bias" = "",
      "T" = "U",
      "_p" = " \\(p",
      "1$" = "1\\)",
      "2$" = "2\\)",
      "3$" = "3\\)",
      "aves" = "bird",
      "camel" = "camelid",
      "yangbat" = "yangochiroptera",
      "yinbat" = "yinpterochiroptera"
    )))

  plot_df_summ <- plot_df %>%
    melt(id.vars = names(plot_df)[1]) %>%
    group_by(!!sym(names(plot_df)[1]), variable) %>%
    summarise(
      median = median(value),
      lower = quantile(value, probs = .025),
      upper = quantile(value, probs = .975)
    )

  assign(
    paste0("SF6_spike_spike_", i),
    ggplot(
      plot_df_summ,
      aes(x = !!sym(names(plot_df)[1]))
    ) +
      geom_ribbon(aes(fill = variable, ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_manual(values = cbbPalette_ordered_bw, name = "host") +
      scale_fill_manual(values = cbbPalette_ordered_bw, name = "host") +
      ylab("Probability") +
      guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
      theme_bw()
  )
}

for (i in 1:length(vars_to_pd_wg)) {
  plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_wg[i]]]) %>%
    bind_rows() %>%
    rename_with(., ~ str_replace_all(., c(
      "_Bias" = "",
      "T" = "U",
      "_p" = " \\(p",
      "1$" = "1\\)",
      "2$" = "2\\)",
      "3$" = "3\\)",
      "aves" = "bird",
      "camel" = "camelid",
      "yangbat" = "yangochiroptera",
      "yinbat" = "yinpterochiroptera"
    )))

  plot_df_summ <- plot_df %>%
    melt(id.vars = names(plot_df)[1]) %>%
    group_by(!!sym(names(plot_df)[1]), variable) %>%
    summarise(
      median = median(value),
      lower = quantile(value, probs = .025),
      upper = quantile(value, probs = .975)
    )

  assign(
    paste0("SF7_spike_wg_", i),
    ggplot(
      plot_df_summ,
      aes(x = !!sym(names(plot_df)[1]))
    ) +
      geom_ribbon(aes(fill = variable, ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_manual(values = cbbPalette_ordered_bw, name = "host") +
      scale_fill_manual(values = cbbPalette_ordered_bw, name = "host") +
      ylab("Probability") +
      guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
      theme_bw()
  )
}

# Supplementary Data Files

# Supplementary Data File 1, all spike protein sequences and their predictions
validate_used %>%
  bind_rows() %>%
  add_column(prediction = predict_class_test) %>%
  cbind(round(predict_prob_test, 5)) %>%
  left_join(df_tax_matcher, by = "accessionversion") %>%
  select(genus, childtaxa_name, accessionversion, name_txt, outcome, prediction, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, metadata_host = name_txt, host_category = outcome, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  mutate_at(vars(host_category, prediction), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  ))) %>%
  write.csv("figures_tables\\Supp_Data_1.csv")

# Supplementary Data File 3, zoonotic spike protein sequences and their predictions
validate_zoonotic %>%
  cbind(round(predict_prob_zoonotic %>% select(aves:yinbat), 5)) %>%
  select(genus, childtaxa_name, accessionversion, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  arrange(virus) %>%
  write.csv("figures_tables\\Supp_Data_3.csv")

# Statistical tests

# One-sample t-tests against null for RSCU values
spike_au_ending <- model_df %>%
  as.data.frame() %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][A|U]$")) %>%
  janitor::remove_constant() %>%
  select(-UAA, -UGA) %>%
  map_df(~ broom::tidy(t.test(., mu = 1, alternative = "greater")), .id = "bias") %>%
  mutate(adj.p = p.adjust(p.value, method = "bonferroni")) %>%
  data.frame() %>%
  arrange(adj.p) %>%
  write.csv("figures_tables\\onesample_t_spike_au_ending.csv")
spike_gc_ending <- model_df %>%
  as.data.frame() %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][G|C]$")) %>%
  janitor::remove_constant() %>%
  select(-UAG) %>%
  map_df(~ broom::tidy(t.test(., mu = 1, alternative = "less")), .id = "bias") %>%
  mutate(adj.p = p.adjust(p.value, method = "bonferroni")) %>%
  data.frame() %>%
  arrange(adj.p) %>%
  write.csv("figures_tables\\onesample_t_spike_gc_ending.csv")

#############################
# Load WGS analysis results #
#############################

load("outputs\\listresults_ml_wg_02_10_20.RData")

# Set unclassified genus
model_df %<>% replace_na(list(genus = "unclassified"))

# Extract parameter optimisation
gridsearch <- lapply(rf_list, function(x) {
  x$results
}) %>%
  bind_rows() %>%
  mutate(
    min.node.size = factor(min.node.size, levels = c("5", "12.5", "20")),
    mtry = factor(mtry, levels = c("5", "12.5", "20"))
  )

# Select holdout used
validate_used <- lapply(unique(model_df$taxid), function(x) {
  model_df %>% filter(taxid == x)
})

# Extract class predictions (excluding zoonotic sequences)
predict_class_test <- Map(function(model, newdata) {
  if (nrow(newdata) > 0) {
    levels(model_df$outcome)[predict(model, newdata = newdata, type = "raw")]
  }
},
model = rf_list,
newdata = validate_used
) %>%
  unlist() %>%
  as.factor()

# Extract probabilistic predictions (excluding zoonotic sequences)
predict_prob_test <- Map(function(model, newdata) {
  if (nrow(newdata) > 0) {
    predict(model, newdata = newdata, type = "prob")
  }
},
model = rf_list,
newdata = validate_used
) %>% bind_rows()

# Calculate confusion matrix
matrix_test <- confusionMatrix(predict_class_test, validate_used %>% bind_rows() %>% pull(outcome) %>% droplevels())

# Calculate one-vs-all confusion matrices
# Sourced from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/ (Matthias Döring, 2018)
matrix_one_vs_all <- vector("list", length(levels(model_df$outcome)))
for (i in seq_along(matrix_one_vs_all)) {
  positive.class <- levels(model_df$outcome)[i]
  # in the i-th iteration, use the i-th class as the positive class
  matrix_one_vs_all[[i]] <- confusionMatrix(predict_class_test, validate_used %>% bind_rows() %>% pull(outcome) %>% droplevels(),
    positive = positive.class
  )
}


# Collect equivalent test set for zoonotic sequences (i.e. MERS-CoV, SARS-CoV, SARS-CoV-2) sampled from humans
set.seed(1223)
validate_zoonotic <- model_df_predownsample %>%
  filter(taxid %in% c(1335626, 694009, 627442, 228404, 2697049, 1263720, 1235996, 1298362, 1306931, 627439, 166124) & outcome == "human") %>%
  split(.$taxid) %>%
  lapply(function(x) {
    x %>% sample_n(ifelse(nrow(x) < 20, nrow(x), 20))
  }) %>%
  bind_rows() %>%
  mutate(outcome = factor(outcome, levels = levels(model_df$outcome)))

# Extract probabilistic predictions for zoonotic sequences (grand mean probabilities over all random forests)
predict_prob_zoonotic <- lapply(rf_list, function(model) {
  validate_zoonotic %>%
    select(taxid, accessionversion) %>%
    cbind(predict(model, newdata = validate_zoonotic, type = "prob"))
}) %>%
  bind_rows() %>%
  group_by(taxid, accessionversion) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(flag = case_when(
    taxid %in% c(1335626, 1263720, 1235996, 1298362, 1306931) ~ "MERS-CoV",
    taxid == 2697049 ~ "SARS-CoV-2",
    taxid %in% c(694009, 627442, 228404) ~ "SARS-CoV"
  ))


# Calculate variable importance
varimp_wg <- lapply(rf_list, function(x) {
  varImp(x)$importance %>%
    rownames_to_column("name") %>%
    mutate(Overall = Overall / 100) %>%
    rename(relGini = Overall)
}) %>%
  bind_rows() %>%
  mutate(name = str_replace_all(name, c(
    "_Bias" = "",
    "T" = "U",
    "_p" = " \\(p",
    "1$" = "1\\)",
    "2$" = "2\\)",
    "3$" = "3\\)"
  )))

varimp_order_wg <- varimp_wg %>%
  group_by(name) %>%
  summarise(mean = mean(relGini)) %>%
  mutate(name = fct_reorder(name, -mean)) %>%
  pull(name) %>%
  levels()

# Read in partial dependence profiles
load("outputs\\listresults_pd_wg_03_10_20.RData")

################################
# Save/plot tables and figures #
################################

# Tables

# Table 1, descriptive stats - mean, SD over each coronavirus species and unranked subspecies
enc_per_spp <- model_df_predownsample %>%
  filter(accessionversion %in% model_df$accessionversion) %>%
  group_by(genus, accessionversion) %>%
  summarise(spp_mean = mean(enc)) %>%
  ungroup()

bind_rows(
  enc_per_spp %>%
    group_by(genus) %>%
    summarise(
      mean = mean(spp_mean),
      sd = sd(spp_mean)
    ),
  data.frame(genus = "Total", mean = mean(enc_per_spp$spp_mean), sd = sd(enc_per_spp$spp_mean))
) %>%
  write.csv("figures_tables\\Table_1_wg.csv")

# Table 2, predictive performance of random forest models
matrix_test$overall %>%
  round(., 3) %>%
  t() %>%
  cbind(.,
    AUC = multiclass.roc(response = validate_used %>%
      bind_rows() %>%
      pull(outcome), predictor = predict_prob_test) %>%
      .$auc %>%
      as.numeric() %>%
      round(3),
    micro_F1 = matrix_one_vs_all %>% get.micro.f1() %>% round(3),
    macro_F1 = matrix_one_vs_all %>% get.macro.f1() %>% round(3)
  ) %>%
  write.csv("figures_tables\\Table_2_wg.csv")

# Supplementary Tables

# Supplementary Table S2, number of sequences per host category
model_df %>%
  group_by(outcome) %>%
  summarise(n_sequences = n(), n_species = n_distinct(taxid)) %>%
  write.csv("figures_tables\\Supp_Table_2_wg.csv")

# Supplementary Table 7, precision/recall by class
matrix_test$byClass %>% write.csv("figures_tables\\Supp_Table_7.csv")

# Figures

# Figure 1, heatmap
png("figures_tables\\Figure_1B.png", width = 12, height = 9, units = "in", res = 320)
heatmap_wgs <- model_df %>%
  as.data.frame() %>%
  set_rownames(.$accessionversion) %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][A|C|G|U]$")) %>%
  janitor::remove_constant() %>%
  as.matrix() %>%
  t() %>%
  heatmap.2(
    density.info = "none",
    trace = "none",
    margins = c(1, 5),
    dendrogram = "col",
    Rowv = "NA",
    labCol = NA,
    ColSideColors = model_df %>% mutate(colsidecol = case_when( # Set side colours using same genus colours as ggplots elsewhere
      genus == "Alphacoronavirus" ~ "#D55E00",
      genus == "Betacoronavirus" ~ "#F0E442",
      genus == "Deltacoronavirus" ~ "#0072B2",
      genus == "Gammacoronavirus" ~ "#009E73",
      genus == "unclassified" ~ "#999999",
    )) %>% pull(colsidecol),
    cexRow = 1.1,
    col = colorRampPalette(c("dodgerblue", "gray10", "firebrick2"))(n = 65),
    breaks = c(seq(0, 0.95, length = 25), seq(0.951, 1.05, length = 6), seq(1.051, 3.7, length = 35)),
    keysize = 0.75, key.par = list(cex = 0.5), key.title = NA, key.xlab = "RSCU",
    lhei = c(2, 10), lwid = c(2, 10)
  )
par(lend = 1) # square line ends for the color legend
legend(
  cex = 0.8, x = -0.05, y = 0.9, xpd = TRUE,
  legend = model_df$genus %>% unique() %>% sort(), fill = c("#D55E00", "#F0E442", "#0072B2", "#009E73", "#999999"), ncol = 1
)
dev.off()

# Compare dendogram clustering to genus using Normalised Mutal Information
heatmap_wgs %>%
  .$colDendrogram %>%
  as.hclust() %>%
  cutree(k = 10) %>%
  as.data.frame() %>%
  set_colnames("cluster") %>%
  rownames_to_column("accessionversion") %>%
  left_join(model_df %>% select(accessionversion, genus)) %>%
  filter(genus != "unclassified") %>%
  with(., NMI(cluster, genus))

# Figure 2, model predictions for held-out coronavirus sequences
valid_df_raw <- data.frame(
  host_category = validate_used %>% bind_rows() %>% pull(outcome),
  accessionversion = validate_used %>% bind_rows() %>% pull(accessionversion),
  genus = validate_used %>% bind_rows() %>% pull(genus),
  childtaxa_name = validate_used %>% bind_rows() %>% pull(childtaxa_name),
  predict_prob_test
)

valid_df_ordered_wg <- rbind(
  valid_df_raw %>% filter(host_category == "aves") %>% arrange(-aves),
  valid_df_raw %>% filter(host_category == "camel") %>% arrange(-camel),
  valid_df_raw %>% filter(host_category == "carnivore") %>% arrange(-carnivore),
  valid_df_raw %>% filter(host_category == "human") %>% arrange(factor(childtaxa_name, levels = valid_df_raw %>%
    filter(host_category == "human") %>%
    arrange(-human) %>%
    pull(childtaxa_name) %>%
    unique()), -human),
  valid_df_raw %>% filter(host_category == "rodent") %>% arrange(-rodent),
  valid_df_raw %>% filter(host_category == "swine") %>% arrange(-swine),
  valid_df_raw %>% filter(host_category == "yangbat") %>% arrange(-yangbat),
  valid_df_raw %>% filter(host_category == "yinbat") %>% arrange(-yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  select(-genus) %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))

F2_wg <- valid_df_ordered_wg %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus HKU1") %>% pull(index) %>% mean(), y = 0.05, label = "HKU1", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus NL63") %>% pull(index) %>% mean(), y = 0.05, label = "NL63", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus 229E") %>% pull(index) %>% mean(), y = 0.05, label = "229E", fill = NULL)
  ) +
  geom_rect(
    data = data.frame(host_category = "human"),
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% min() - 0.5,
      xmax = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% max() + 0.5,
      ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = data.frame(host_category = "human"), size = 3,
    aes(x = valid_df_ordered_wg %>% filter(childtaxa_name == "Human coronavirus OC43") %>% pull(index) %>% mean(), y = 0.05, label = "OC43", fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Figure 3, model predictions for zoonotic coronavirus sequences
F3_wg <- rbind(
  predict_prob_zoonotic %>% filter(flag == "MERS-CoV") %>% arrange(-camel),
  predict_prob_zoonotic %>% filter(flag == "SARS-CoV") %>% arrange(-rodent),
  predict_prob_zoonotic %>% filter(flag == "SARS-CoV-2") %>% arrange(-yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  melt(id.vars = c("accessionversion", "taxid", "flag")) %>%
  mutate_at(vars(variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  ))) %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_grid(~flag, scales = "free_x", space = "free_x")

# Figure 4, variable importance
varimp_comp <- varimp_wg %>%
  group_by(name) %>%
  summarise(mean_wgs = mean(relGini)) %>%
  left_join(varimp_spike %>%
    group_by(name) %>%
    summarise(mean_spike = mean(relGini)),
  by = "name"
  ) %>%
  mutate(type = case_when(
    grepl("^[A|C|G|U]$", name) ~ "nucleotide biases",
    grepl("^[A|C|G|U][A|C|G|U] ", name) ~ "dinucleotide biases",
    grepl("^[A|C|G|U][A|C|G|U][A|C|G|U]$", name) ~ "codon biases"
  ))

F4 <- ggplot(
  varimp_comp,
  aes(x = mean_wgs, y = mean_spike, colour = type)
) +
  geom_point(size = 4, alpha = 0.5) +
  scale_colour_manual(values = c("#F0E442", "#CC79A7", "#56B4E9"), name = "Predictor type") +
  xlab("Mean rel. Gini decrease, whole genome") +
  ylab("Mean rel. Gini decrease, spike") +
  ggrepel::geom_text_repel(
    data = varimp_comp %>% filter(name %in% c(varimp_comp %>% slice_max(mean_wgs, n = 10) %>% pull(name), varimp_comp %>% slice_max(mean_spike, n = 10) %>% pull(name))),
    aes(label = name), colour = "black"
  ) +
  theme_bw() +
  theme(legend.justification = c(1, 1), legend.position = c(.98, .98), panel.spacing = unit(1.1, "lines"), legend.title = element_blank())

F4 + ggsave("figures_tables\\Figure_4.png", width = 6, height = 6)

# Supplementary Figures

# Supplementary Figure S2 - grid search results
SF2_wg <- ggplot(gridsearch, aes(y = Accuracy, x = min.node.size, fill = mtry)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Cross-validation (inner loop) accuracy")

# Supplementary Figure S3, dinucleotide bias values by position
SF3_wg <- model_df %>%
  rename_with(., ~ str_replace_all(., c(
    "_Bias" = "",
    "T" = "U",
    "_p" = " \\(p",
    "1$" = "1\\)",
    "2$" = "2\\)",
    "3$" = "3\\)"
  ))) %>%
  select(matches("^[A|C|G|U][A|C|G|U] ")) %>%
  melt() %>%
  mutate(
    dinuc = substr(variable, 1, 2),
    position = substr(variable, 6, 6),
    position = case_when(
      position == "1" ~ "1-2",
      position == "2" ~ "2-3",
      position == "3" ~ "3-1"
    )
  ) %>%
  ggplot(aes(y = value, x = position, fill = position)) +
  geom_boxplot(outlier.shape = NA, width = 0.9) +
  facet_wrap(. ~ variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  ylab("Dinucleotide bias") +
  geom_hline(yintercept = 1, alpha = 0.3, color = "grey50", lty = "dashed", size = 1) +
  facet_wrap(~dinuc)

# Supplementary Figure S4, model predictions for held-out coronavirus sequences blocked by genus
valid_df_ordered_wg_genus <- rbind(
  valid_df_raw %>% filter(host_category == "aves") %>% arrange(genus, -aves),
  valid_df_raw %>% filter(host_category == "camel") %>% arrange(genus, -camel),
  valid_df_raw %>% filter(host_category == "carnivore") %>% arrange(genus, -carnivore),
  valid_df_raw %>% filter(host_category == "human") %>% arrange(genus, -human),
  valid_df_raw %>% filter(host_category == "rodent") %>% arrange(genus, -rodent),
  valid_df_raw %>% filter(host_category == "swine") %>% arrange(genus, -swine),
  valid_df_raw %>% filter(host_category == "yangbat") %>% arrange(genus, -yangbat),
  valid_df_raw %>% filter(host_category == "yinbat") %>% arrange(genus, -yinbat)
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index", "genus")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))


known_combs <- valid_df_ordered_wg_genus %>% distinct(host_category, genus)
wg_genus_annots <- list()
for (i in 1:nrow(known_combs)) {
  wg_genus_annots[[i]] <- data.frame(
    host_category = known_combs$host_category[i],
    genus = known_combs$genus[i],
    xmin = valid_df_ordered_wg_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% min() - 0.5,
    xmax = valid_df_ordered_wg_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% max() + 0.5,
    xmean = valid_df_ordered_wg_genus %>%
      filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>%
      pull(index) %>% mean()
  ) %>%
    mutate(label = case_when(
      genus == "Alphacoronavirus" ~ "alpha",
      genus == "Betacoronavirus" ~ "beta",
      genus == "Gammacoronavirus" ~ "gamma",
      genus == "Deltacoronavirus" ~ "delta",
      genus == "unclassified" ~ "U"
    ))
}

wg_genus_annots <- bind_rows(wg_genus_annots)

SF4_wg_genus <- valid_df_ordered_wg_genus %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = wg_genus_annots,
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = xmin, xmax = xmax, ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.5
  ) +
  geom_text(
    data = wg_genus_annots, size = 3,
    aes(x = xmean, y = 0.05, label = label, fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Supplementary Figure S5, model predictions for held-out coronavirus sequences blocked by species
valid_df_ordered_wg_spp <- rbind(
  rearrange_to_plot_pred_fig_spp("aves"),
  rearrange_to_plot_pred_fig_spp("camel"),
  rearrange_to_plot_pred_fig_spp("carnivore"),
  rearrange_to_plot_pred_fig_spp("human"),
  rearrange_to_plot_pred_fig_spp("rodent"),
  rearrange_to_plot_pred_fig_spp("swine"),
  rearrange_to_plot_pred_fig_spp("yangbat"),
  rearrange_to_plot_pred_fig_spp("yinbat")
) %>%
  mutate(accessionversion = factor(accessionversion, levels = accessionversion)) %>%
  group_by(host_category) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  select(-genus) %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  )))


known_combs <- valid_df_ordered_wg_spp %>% distinct(host_category, childtaxa_name)
wg_spp_annots <- list()
for (i in 1:nrow(known_combs)) {
  wg_spp_annots[[i]] <- data.frame(
    host_category = known_combs$host_category[i],
    childtaxa_name = known_combs$childtaxa_name[i],
    xmin = valid_df_ordered_wg_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% min() - 0.5,
    xmax = valid_df_ordered_wg_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% max() + 0.5,
    xmean = valid_df_ordered_wg_spp %>%
      filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>%
      pull(index) %>% mean()
  )
}

wg_spp_annots <- bind_rows(wg_spp_annots) %>%
  left_join(spike_spp_annots %>% select(childtaxa_name, label)) # use labels as in panel A

SF5_wg_spp <- valid_df_ordered_wg_spp %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  geom_rect(
    data = wg_spp_annots,
    aes(
      x = NULL, y = NULL, fill = NULL,
      xmin = xmin, xmax = xmax, ymin = 0, ymax = 1
    ), colour = "black", fill = NA, lwd = 0.3
  ) +
  geom_text(
    data = wg_spp_annots, size = 1.4, angle = 90,
    aes(x = xmean, y = -0.05, label = label, fill = NULL)
  ) +
  scale_y_continuous(expand = c(0, 0.05)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~host_category, scales = "free_x", nrow = 2)

# Supplementary Figure S6 & S7, partial dependence
for (i in 1:length(vars_to_pd_spike)) {
  plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_spike[i]]]) %>%
    bind_rows() %>%
    rename_with(., ~ str_replace_all(., c(
      "_Bias" = "",
      "T" = "U",
      "_p" = " \\(p",
      "1$" = "1\\)",
      "2$" = "2\\)",
      "3$" = "3\\)",
      "aves" = "bird",
      "camel" = "camelid",
      "yangbat" = "yangochiroptera",
      "yinbat" = "yinpterochiroptera"
    )))

  plot_df_summ <- plot_df %>%
    melt(id.vars = names(plot_df)[1]) %>%
    group_by(!!sym(names(plot_df)[1]), variable) %>%
    summarise(
      median = median(value),
      lower = quantile(value, probs = .025),
      upper = quantile(value, probs = .975)
    )

  assign(
    paste0("SF6_wg_spike_", i),
    ggplot(
      plot_df_summ,
      aes(x = !!sym(names(plot_df)[1]))
    ) +
      geom_ribbon(aes(fill = variable, ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_manual(values = cbbPalette_ordered_bw, name = "host") +
      scale_fill_manual(values = cbbPalette_ordered_bw, name = "host") +
      ylab("Probability") +
      guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
      theme_bw()
  )
}

for (i in 1:length(vars_to_pd_wg)) {
  plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_wg[i]]]) %>%
    bind_rows() %>%
    rename_with(., ~ str_replace_all(., c(
      "_Bias" = "",
      "T" = "U",
      "_p" = " \\(p",
      "1$" = "1\\)",
      "2$" = "2\\)",
      "3$" = "3\\)",
      "aves" = "bird",
      "camel" = "camelid",
      "yangbat" = "yangochiroptera",
      "yinbat" = "yinpterochiroptera"
    )))

  plot_df_summ <- plot_df %>%
    melt(id.vars = names(plot_df)[1]) %>%
    group_by(!!sym(names(plot_df)[1]), variable) %>%
    summarise(
      median = median(value),
      lower = quantile(value, probs = .025),
      upper = quantile(value, probs = .975)
    )

  assign(
    paste0("SF7_wg_wg_", i),
    ggplot(
      plot_df_summ,
      aes(x = !!sym(names(plot_df)[1]))
    ) +
      geom_ribbon(aes(fill = variable, ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_manual(values = cbbPalette_ordered_bw, name = "host") +
      scale_fill_manual(values = cbbPalette_ordered_bw, name = "host") +
      ylab("Probability") +
      guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
      theme_bw()
  )
}

# Supplementary Data Files

# Supplementary Data File 2, all whole genome sequences and their predictions
validate_used %>%
  bind_rows() %>%
  add_column(prediction = predict_class_test) %>%
  cbind(round(predict_prob_test, 5)) %>%
  left_join(df_tax_matcher, by = "accessionversion") %>%
  select(genus, childtaxa_name, accessionversion, name_txt, outcome, prediction, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, metadata_host = name_txt, host_category = outcome, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  mutate_at(vars(host_category, prediction), funs(str_replace_all(
    .,
    c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")
  ))) %>%
  write.csv("figures_tables\\Supp_Data_2.csv")

# Supplementary Data File 4, zoonotic whole genome sequences and their predictions
validate_zoonotic %>%
  cbind(round(predict_prob_zoonotic %>% select(aves:yinbat), 5)) %>%
  select(genus, childtaxa_name, accessionversion, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  arrange(virus) %>%
  write.csv("figures_tables\\Supp_Data_4.csv")

# Statistical tests

# One-sample t-tests against null for RSCU values
wgs_au_ending <- model_df %>%
  as.data.frame() %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][A|U]$")) %>%
  janitor::remove_constant() %>%
  select(-UAA, -UGA) %>%
  map_df(~ broom::tidy(t.test(., mu = 1, alternative = "greater")), .id = "bias") %>%
  mutate(adj.p = p.adjust(p.value, method = "bonferroni")) %>%
  data.frame() %>%
  arrange(adj.p) %>%
  write.csv("figures_tables\\onesample_t_wgs_au_ending.csv")
wgs_gc_ending <- model_df %>%
  as.data.frame() %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][G|C]$")) %>%
  janitor::remove_constant() %>%
  select(-UAG) %>%
  map_df(~ broom::tidy(t.test(., mu = 1, alternative = "less")), .id = "bias") %>%
  mutate(adj.p = p.adjust(p.value, method = "bonferroni")) %>%
  data.frame() %>%
  arrange(adj.p) %>%
  write.csv("figures_tables\\onesample_t_wgs_gc_ending.csv")

##########################################################
# Combined figure panels for spike protein, wgs analysis #
##########################################################

F2 <- F2_spike + (F2_wg + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2) &
  theme(legend.position = "right")

F2 + ggsave("figures_tables\\Figure_2.png", width = 15, height = 9)

F3 <- F3_spike + (F3_wg + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2) &
  theme(legend.position = "right")

F3 + ggsave("figures_tables\\Figure_3.png", width = 10.5, height = 5)

SF2 <- SF2_spike + (SF2_wg + guides(fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

SF2 + ggsave("figures_tables\\Supp_Figure_S2.png", width = 10, height = 4)

SF3 <- SF3_spike + (SF3_wg + guides(fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

SF3 + ggsave("figures_tables\\Supp_Figure_S3.png", width = 10, height = 6)

SF4 <- SF4_spike_genus + (SF4_wg_genus + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2) &
  theme(legend.position = "right")

SF4 + ggsave("figures_tables\\Supp_Figure_S4.png", width = 15, height = 9)

SF5 <- SF5_spike_spp + (SF5_wg_spp + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2) &
  theme(legend.position = "right")

SF5 + ggsave("figures_tables\\Supp_Figure_S5.png", width = 17, height = 9)

spike_spp_annots %>%
  select(label, childtaxa_name) %>%
  distinct() %>%
  rename(virus = childtaxa_name) %>%
  write.table("figures_tables\\Supp_Figure_S5_labels.txt", sep = "   ", row.names = FALSE)

SF6 <- SF6_spike_spike_1 +
  (SF6_spike_spike_2 + guides(colour = "none", fill = "none")) +
  (SF6_spike_spike_3 + guides(colour = "none", fill = "none")) +
  (SF6_spike_spike_4 + guides(colour = "none", fill = "none")) +
  (SF6_wg_spike_1 + guides(colour = "none", fill = "none")) +
  (SF6_wg_spike_2 + guides(colour = "none", fill = "none")) +
  (SF6_wg_spike_3 + guides(colour = "none", fill = "none")) +
  (SF6_wg_spike_4 + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

SF6 + ggsave("figures_tables\\Supp_Figure_S6.png", width = 16, height = 8)

SF7 <- SF7_wg_wg_1 +
  (SF7_wg_wg_2 + guides(colour = "none", fill = "none")) +
  (SF7_wg_wg_3 + guides(colour = "none", fill = "none")) +
  (SF7_wg_wg_4 + guides(colour = "none", fill = "none")) +
  (SF7_spike_wg_1 + guides(colour = "none", fill = "none")) +
  (SF7_spike_wg_2 + guides(colour = "none", fill = "none")) +
  (SF7_spike_wg_3 + guides(colour = "none", fill = "none")) +
  (SF7_spike_wg_4 + guides(colour = "none", fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

SF7 + ggsave("figures_tables\\Supp_Figure_S7.png", width = 16, height = 8)
