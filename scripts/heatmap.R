h_map_df <- eid2_full_df %>% 
  select(family_name, names(eid2_full_df)[4254:4267]) %>%
  mutate_at(vars(-family_name), . %>% as.character() %>% as.numeric()) %>%
  group_by(family_name) %>% 
  summarise_all(funs(sum)) %>%
  melt() %>%
  rename(order = variable, count = value)

# Calculate number of human-shared per virus/host group combination and merge
h_map_df <- merge(h_map_df,  
                  eid2_full_df %>% 
                    select(family_name, names(eid2_full_df)[4254:4267]) %>%
                    mutate_at(vars(-family_name), . %>% as.character() %>% as.numeric()) %>%
                    filter(human == 1) %>%
                    group_by(family_name, .drop=FALSE) %>% 
                    summarise_all(funs(sum)) %>%
                    melt() %>%
                    rename(order = variable, num = value),
                    by = c("family_name", "order"), all.x=TRUE) %>%
  mutate(num = replace_na(num, 0)) %>%
  mutate(prop = num/count) %>%
  mutate(plot = ifelse(order == "human", "x_human", "mammal")) %>% # Plot humans separately
  mutate(order = fct_relevel(order, rev)) %>% # Reorder levels for plot
  mutate(family_name = vir_family_short(family_name))

ggplot(h_map_df, aes(family_name, order, fill = count)) + 
  geom_tile(colour = "gray80", size=0.5) +  
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_discrete(expand = c(0, 0.5)) +
  scale_fill_gradient(low = "white", high = "firebrick3", name="Virus count", limits = c(0,25)) +
  theme(panel.margin.y=unit(1, "lines"), axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y = element_text(size=12), axis.title = element_text(size = 14), legend.text = element_text(size = 11), legend.title=element_text(size=12), strip.text = element_blank(), strip.background = element_blank()) +
  facet_grid(plot ~ ., scales="free",drop=T,space="free") +
  xlab("Virus family") + 
  ylab("Host group")

# ENC heatmap - TO FINISH
merge(merge(vir_enc_df, 
      accessions_table %>% select(genaccessionref, taxid) %>% unique,
      by.x = "SeqName",
      by.y = "genaccessionref") %>% select(-SeqName),
      eid2_full_df,
      by.x = "taxid",
      by.y = "chosen_refseq_taxid") %>%
 select(family_name, names(eid2_full_df)[4254:4267], mean_enc) %>%
  melt(id=family) %>%
  
  
  
  mutate_at(vars(-family_name), . %>% as.character() %>% as.numeric()) %>%
  group_by() %>% 
  summarise_all(funs(sum)) %>%
  melt() %>%
  rename(order = variable, count = value)

#Kind of working but some values wildly in the wrong place

ggplot(h_map_df %>% filter(order != "human") %>% arrange(family_name, order), aes(family_name, order, fill = prop)) + 
  geom_tile(colour="gray80", size=0.5) + 
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_discrete(expand = c(0, 0.5)) +
  geom_text(data=h_map_df %>% filter(order != "human" & count > 0), 
            aes(label=num, x=as.numeric(as.factor(family_name))-0.25, y=as.numeric(as.factor(order))+0.25, size=14)) +
  geom_text(data=h_map_df %>% filter(order != "human" & count > 0), 
            aes(label="of", size=14)) +
  geom_text(data=h_map_df %>% filter(order != "human" & count > 0), 
            aes(label=count, x=as.numeric(as.factor(family_name))+0.25, y=as.numeric(as.factor(order))-0.25, size=14)) +
  scale_fill_gradient(low = "white",  high = "red", name="prop. zoon", na.value="gray20") +
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=12), axis.text.y = element_text(size=12), axis.title = element_text(size = 14), legend.text = element_text(size = 11), legend.title=element_text(size=12), strip.text = element_blank(), strip.background = element_blank())  +
#  facet_grid(plot ~ ., scales="free",drop=T,space="free") +
  xlab("Virus family") + 
  ylab("Host group")