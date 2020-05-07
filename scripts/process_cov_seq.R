###############################
# Extract sequences from NCBI #
###############################

# Extract all viruses in family Coronaviridae and their taxonomic IDs
allcov_df <- downstream("Coronaviridae", db="ncbi")$Coronaviridae %>% filter(rank %in% c("no rank", "species"))

# Find and obtain results based on taxonomy IDs of all coronaviruses
if (load_prev_seqs == TRUE){
  load(file="data\\allcov_results.RData")
} else {
  allcov_results <- lapply(as.list(allcov_df$childtaxa_id), function(x) Seq_searcher(x))
  save(allcov_results, file="data\\allcov_results.RData")
}

# Load IDs of sequences to exclude as they're known a priori to have no CDS annotation
load("data//missingaccession.RData")

# Append genus name, sequence counts of viruses found
allcov_df <- data.frame(allcov_df, 
                        n_seqs = lapply(allcov_results, function(x) x$ids %>% subset(!(. %in% (missinguid %>% as.character))) %>% length) %>% unlist())

genus_list <- lapply(allcov_df$childtaxa_id %>% classification, function(x) x %>% filter(rank == "genus") %>% pull(name))
is.na(genus_list) <- lengths(genus_list) == 0

allcov_df$genus <- genus_list %>% unlist %>% unname

# Extract IDs of individual sequences
allcov_seqids <- lapply(allcov_results, function(x) x$ids %>% subset(!(. %in% (missinguid %>% as.character)))) %>% Filter(length, .)

# Fetch metadata information for those sequences (title, accession). NB entrez_summary doesn't accept empty objects, so filtering out those IDs with no sequences
if (load_prev_seqs == TRUE){
  load(file="data\\allcov_found_summaries.RData")
} else {
  allcov_found_summaries <- pblapply(allcov_seqids, function(x) Seq_summary(x))
  save(allcov_found_summaries, file="data\\allcov_found_summaries.RData")
}

# Fetch sequences as FASTA
if (load_prev_seqs == TRUE){
  load(file="data\\allcov_FASTA.RData")
} else {
  allcov_FASTA <- pblapply(allcov_seqids, function(x) Seq_FASTA(x))
  save(allcov_FASTA, file="data\\allcov_FASTA.RData")
}

# Write FASTA file
write(unlist(allcov_FASTA), file="data\\allcov_GenBank.fasta")

# Identify n sequences for SARS-CoV-2
allcov_df %>% filter(childtaxa_id == 2697049)

# NB numbers don't match between this and SARS-CoV-2 special page on GenBank because some of them have incorrectly just been filed under taxid 694009 for SARS-CoV parent species...
# Load in FASTA from manually downloaded GISAID repository (2020-03-27 18:15:04 GMT)
sarscov2_GenBank_FASTA <- read.fasta("data\\sarscov2_GenBank.fasta")
sarscov2_GISAID_FASTA <- read.fasta("data\\sarscov2_GISAID.fasta")

length(sarscov2_GenBank_FASTA) # total CDS
n_man_genbank_seqs <- sarscov2_GenBank_FASTA[lapply(sarscov2_GenBank_FASTA, function(x) x %>% 
                                                      attributes %>% 
                                                      .$Annot %>% 
                                                      grepl("gene=S", .)) %>% unlist] %>% length # total S sequences - should be larger but it's smaller?? some might be missing annotation?
length(sarscov2_GISAID_FASTA) # total complete genomes

########################
# Add host information #
########################

eid2_cov <- read.csv("data\\species_species_dataset_aves_mammals.csv", fileEncoding="UTF-8-BOM")  %>% 
  filter(PathogenTaxID %in% allcov_df$childtaxa_id) %>%
  distinct(PathogenTaxID, HostTaxID, .keep_all = TRUE) # remove duplicated host-virus pairs

# Recast data to create 0/1 outcome variables as to whether each virus is known to infect each host type
eid2_outcomes <- dcast(eid2_cov, PathogenTaxID ~Hc1, length) %>% rename_at(vars(-PathogenTaxID), ~ paste0("h_",gsub(" ","_",.)))
eid2_outcomes %<>% mutate_at(vars(-PathogenTaxID), ~1 * (. > 0))
eid2_outcomes %<>% mutate(eid2_n_groups = rowSums(eid2_outcomes[,-1]))
eid2_outcomes %<>% mutate_at(vars(-PathogenTaxID, -eid2_n_groups), funs(factor))

# Count number of host species infected
# Including the following as host species based on EID2 data structure: humans, domestic subspecies (e.g canis lupus separate to canis lupus familiaris)
eid2_outcomes <- merge(eid2_outcomes,
                       eid2_cov %>% 
                         group_by(PathogenTaxID) %>% 
                         summarise(eid2_n_species = n_distinct(HostTaxID)),
                       by = "PathogenTaxID")

allcov_df <- merge(allcov_df, eid2_outcomes,
                   by.x = "childtaxa_id",
                   by.y = "PathogenTaxID",
                   all.x = TRUE) %>% replace(., is.na(.), 0)

# Specify SARS-CoV-2 as a human virus, which currently isn't recorded in EID2
allcov_df[allcov_df$childtaxa_id == 2697049,"h_human"] <- "1"
allcov_df[allcov_df$childtaxa_id == 2697049,"eid2_n_groups"] <- 1
allcov_df[allcov_df$childtaxa_id == 2697049,"eid2_n_species"] <- 1

##############################
# Clean and process metadata #
##############################

# Extract metadata into df format
allcov_found_summaries %<>% flatten
class(allcov_found_summaries) <-  c("list", "esummary_list") # coerce to esummary list
allcov_meta_df <- sapply(c("uid", "caption", "title", "extra", "gi", "createdate", "updatedate", "flags", "taxid", "slen", "biomol", "moltype", "topology", "sourcedb", "segsetsize", "projectid", "genome", "subtype", "subname", "assemblygi", "assemblyacc", "tech", "completeness", "geneticcode", "strand", "organism", "strain", "biosample", "accessionversion")
                         , function(x) extract_from_esummary(allcov_found_summaries, x) %>% as.vector) %>% as.data.frame
allcov_meta_df$meta_length <- extract_from_esummary(allcov_found_summaries, "statistics") %>% t %>% .[,2] %>% lapply(., function(x) x[1]) %>% unlist

# Clean up metadata and assign each entry as complete or partial genome
allcov_meta_df %<>% mutate(complete = case_when(
  grepl("complete genome", title, ignore.case=TRUE) ~ "whole_genome",    # Assign as whole genome if labelled so
  grepl("partial and complete cds", title, ignore.case=TRUE) ~ "partial_spike",    # Else assign partial if mixed
  metadata_title_cleaner(title) == "partial" ~ "partial_spike",          # Else apply the metadata cleaner
  metadata_title_cleaner(title) == "complete"  ~ "complete_spike",
  grepl("partial", title, ignore.case=TRUE) ~ "partial_spike",           # If the complex metadata cleaner can't assign it (which is only 3% of records) then just do a simple search match
  grepl("complete", title, ignore.case=TRUE) ~ "complete_spike")) %>% 
  replace_na(list(complete = "partial_spike"))                           # Else no mention of completeness (which is only 0.3% of records), then just assume not complete

## Code to assist choosing regexp for assigning complete
# table(allcov_meta_df$completeness, exclude = NULL) #MOST NOT LABELLED
# allcov_meta_df %>% filter(complete == "partial_spike" & completeness == "complete") %>% select(title) # check against metadata values where they exist
# allcov_meta_df %>% filter(complete == "complete_spike" & completeness == "partial") %>% select(title) # looks like the metadata "completeness" field is wrong or at least conflicts with title

# Re-format sub-metadata and add in as columns
sub_meta_df <- apply(allcov_meta_df, 1, function(x) {
  x["subname"] %>% 
    str_split("\\|") %>% 
    unlist %>% 
    t %>% 
    data.frame %>% 
    set_colnames(x["subtype"] %>% 
                   str_split("\\|") %>% 
                   unlist)
}) %>% bind_rows

sub_meta_df <- sub_meta_df[,which(colnames(sub_meta_df) != "")] # Remove blank columns created by empty metadata
sub_meta_df %<>% rename_all(function(x){paste0("meta_", x)}) # Rename to label them as sub metadata
allcov_meta_df %<>% cbind(sub_meta_df)
rm(sub_meta_df)

if (load_prev_seqs == TRUE){
  
  load(file="data\\tax_matcher.RData")
  
} else {
  
  # Setup empty list
  list_tax_matcher <- vector("list", nrow(allcov_meta_df)) 
  
  # Resolve host metadata to taxonomic species
  for (i in 1:nrow(allcov_meta_df)){
    
    list_tax_matcher[[i]] <- 
      if (!(allcov_meta_df %>% slice(i) %>% pull(meta_host) %>% is.na) == T){
        allcov_meta_df %>% 
          slice(i) %>% 
          pull(meta_host) %>% 
          gsub(" \\(.*$|\\;.*$|\\,.*$","", .) %>%
          gsub(" sp\\.| spp\\.","", .) %>% 
          name2taxid(out_type = "summary") %>%         # If host metadata is present, resolve and add in accessionversion
          mutate(accessionversion = allcov_meta_df %>% 
                   slice(i) %>% 
                   pull(accessionversion))
      } else {
        tibble(name_txt = character(), 
               tax_id = character(), 
               accessionversion = factor()) # If host metadata is not present, return an empty tibble so list can still be collapsed with rbindlist
      }
    
    if (list_tax_matcher[[i]] %>% nrow == 0) { # If empty (either no host metadata, or host metadata not found), return given term alongside NA, and accessionversion
      list_tax_matcher[[i]] <- tibble(name_txt = allcov_meta_df %>% 
                                        slice(i) %>% 
                                        pull(meta_host), 
                                      tax_id = NA, 
                                      accessionversion = allcov_meta_df %>% 
                                        slice(i) %>% 
                                        pull(accessionversion))
    }
  }
  
  # Collapse to data frame
  df_tax_matcher <- data.table::rbindlist(list_tax_matcher) 
  
  # Manually resolve remaining unmatched host labels to lowest possible unambiguous taxonomic level
  df_tax_matcher %<>% mutate(
    tax_id = case_when(
      name_txt == "Camel"|name_txt == "camel" ~ "9836",  # Fix common names that match multiple taxonomic ids with their lowest common parent taxid
      name_txt == "mink" ~ "169418",
      name_txt == "rat" ~ "10114",
      name_txt == "mouse" ~ "10088",
      !is.na(tax_id) ~ tax_id,
      grepl("Gallus gallus|poultry|chicken|chikcen|breeder|layer|broiler|\\<hen\\>|\\<hens\\>|chick with 18 days|embryonated egg", name_txt, ignore.case = TRUE) ~ "9031",    # Fix common domestic species
      grepl("Bos taurus|cow|calf", name_txt, ignore.case = TRUE) ~ "9913",
      grepl("Sus scrofa|pig|porcine", name_txt, ignore.case = TRUE) ~ "9823",
      grepl("Equus caballus|horse", name_txt, ignore.case = TRUE) ~ "9706",
      grepl("Felis catus|feline", name_txt, ignore.case = TRUE) ~ "9685",
      grepl("Canis lupus familiaris|canine|dog", name_txt, ignore.case = TRUE) ~ "9615",
      grepl("Equus africanus asinus|donkey", name_txt, ignore.case = TRUE) ~ "9793",
      grepl("Rhinilophus ferrumequinum", name_txt, ignore.case = TRUE) ~ "59479",  # Fix clear instances of errors or undetectable otherwise unambiguous species
      grepl("Chaerephon plicata", name_txt, ignore.case = TRUE) ~ "478698",
      grepl("Pipistrellus cf. hesperidus", name_txt, ignore.case = TRUE) ~ "294653",
      grepl("bottlenose dolphin", name_txt, ignore.case = TRUE) ~ "9738",
      grepl("pangolin", name_txt, ignore.case = TRUE) ~ "9972",
      grepl("guinea fowl", name_txt, ignore.case = TRUE) ~ "8990",
      grepl("sparrow", name_txt, ignore.case = TRUE) ~ "9126",
      grepl("\\<bat\\>", name_txt, ignore.case = TRUE) ~ "9397",
      grepl("shorebird", name_txt, ignore.case = TRUE) ~ "8906",
      grepl("wigeon", name_txt, ignore.case = TRUE) ~ "1526411",
      grepl("white-eye", name_txt, ignore.case = TRUE) ~ "36297",
      grepl("night-heron", name_txt, ignore.case = TRUE) ~ "8899",
      grepl("magpie-robin", name_txt, ignore.case = TRUE) ~ "125862",
      grepl("teal", name_txt, ignore.case = TRUE) ~ "8835",
      grepl("quail", name_txt, ignore.case = TRUE) ~ "9005",
      grepl("Chinese bulbul", name_txt, ignore.case = TRUE) ~ "125283",
      grepl("antelope", name_txt, ignore.case = TRUE) ~ "9895",
      grepl("civet", name_txt, ignore.case = TRUE) ~ "219112",
      grepl("pheasant|peafowl", name_txt, ignore.case = TRUE) ~ "9072",
      grepl("fox", name_txt, ignore.case = TRUE) ~ "9608",
      TRUE ~ tax_id)) %>% 
    distinct
  
  # For each host taxonomic ID retrieve and arrange into data frame class, order, family, genus, species
  host_taxize_full_class <- df_tax_matcher %>% 
    filter(!is.na(tax_id)) %>% 
    pull(tax_id) %>% 
    unique %>%  
    classification(db="ncbi") %>%
    lapply(function(x) x %>%
             data.frame %>%
             filter(rank == "class" | rank == "order" | rank == "family" | rank == "genus" | rank == "species") %>%
             # Reshape taxonomy tables to give one row per host with family/genus/species*id/name columns
             melt(id="rank") %>%
             unite(rank_var, c("rank", "variable")) %>%
             spread(rank_var, value)) %>%
    bind_rows() %>%
    mutate(tax_id = df_tax_matcher %>% filter(!is.na(tax_id)) %>% pull(tax_id) %>% unique)
  
  # Merge in full taxonomic information
  df_tax_matcher %<>% left_join(host_taxize_full_class, by = "tax_id")
  
  ## Code to assist choosing regexp for assigning host metadata
  ## Calculate % of sequences with host metadata successfully matched
  # (df_tax_matcher %>% filter(!is.na(tax_id)) %>% pull(accessionversion) %>% unique %>% length)*100/(allcov_meta_df %>% filter(!is.na(meta_host)) %>% nrow)
  ## Identify host metadata terms not matched successfully
  #allcov_meta_df %>% filter(!is.na(meta_host) & !(accessionversion %in% df_tax_matcher$accessionversion)) %>% pull(meta_host) %>% unique
  ## Return ambiguous names that cannot be resolved
  #unmatched <- df_tax_matcher %>% filter(is.na(tax_id) & !is.na(name_txt)) %>% pull(name_txt) %>% unique
  # # Identify common names matching multiple taxonomic entries 
  # df_tax_matcher %>% group_by(accessionversion) %>% filter(n() > 1) %>% pull(name_txt) %>% unique
  ## Calc lowest common taxid for common names matching multiple taxonomic entries - slow
  #df_tax_matcher %>% group_by(accessionversion) %>% filter(n() > 1) %>% summarise(taxize::lowest_common(tax_id, db = "ncbi") %>% pull(id))
  
  save(df_tax_matcher, file="data\\tax_matcher.RData")
}

# Append full metadata host label taxonomy to sequence-level data
allcov_meta_df %<>% left_join(df_tax_matcher, by = "accessionversion")

# Count host categories per coronavirus and append to species-level data
allcov_df %<>% left_join(allcov_meta_df %>% 
                           group_by(taxid) %>% 
                           summarise_at(c("class_name","order_name","family_name","genus_name","species_name"), n_distinct, na.rm = TRUE) %>%
                           rename_at(vars(-taxid), ~ paste0("n_", gsub("_name","",.))),
                         by = c("childtaxa_id" = "taxid"))


# Read FASTA back in in coRdon format and calculate codon table
cov_cord <- readSet(file="data\\allcov_GenBank.fasta")

# Construct summary df for individual CDS
cov_enc_df <- data.frame(title = cov_cord %>% names,
                         accessionversion = cov_cord %>% names %>% str_match("lcl\\|(.*?)_cds_") %>% .[,2] %>% as.character,
                         gene = cov_cord %>% names %>% str_match("\\[gene=(.*?)\\]") %>% .[,2],
                         protein = cov_cord %>% names %>% str_match("\\[protein=(.*?)\\]") %>% .[,2],
                         length = cov_cord %>% width)

# Extract and append CDS location information by text extraction from title
cov_enc_df %<>% cbind(cov_enc_df %>% 
                        pull(title) %>% 
                        str_match("\\[location=(.*?)\\]") %>% 
                        .[,2] %>% 
                        gsub("join\\(|\\)|>|<","",.) %>% 
                        gsub(",","..",.) %>% 
                        str_split(., "\\..") %>% 
                        lapply(., function(x) x %>% 
                                 t %>%  
                                 as.data.frame) %>% 
                        bind_rows %>% 
                        mutate_all(., as.numeric) %>%
                        rename_all(., ~ gsub("V","loc_",.)))

# Correct frameshift/slippage causing double accounting of elements of some ORF1ab CDS by specifying new location to start calculating metrics at
cov_enc_df %<>% left_join(
  cov_enc_df %>% filter(                                                                                          # Only apply correction where needed
    accessionversion %in% (cov_enc_df %>% filter(grepl("join",title)) %>% pull(accessionversion)) &               # Filter to sequences containing a CDS that join overlapping components
      accessionversion %in% (allcov_meta_df %>% filter(complete == "whole_genome") %>% pull(accessionversion))    # Filter to sequences that describe a whole genome
  ) %>% group_by(accessionversion, loc_1) %>% arrange(accessionversion, loc_1) %>% filter(n()>1) %>%              # Filter to sequences with two different CDS that both start at the same location
    mutate(start_new = ifelse(is.na(loc_3), lag(loc_2) - (loc_1 - 1) + 1, 1)) %>%                                 # For the shorter sequence of ORF1a, set new start location based on slippage location of ORF1b
    ungroup() %>%
    select(title, start_new),
  by = "title") %>%
  replace_na(list(start_new = 1))

# Apply new CDS start locations
cov_cord %<>% subseq(start = cov_enc_df$start_new)

# Calculate genomic composition metrics for individual CDS
cov_enc_df %<>% cbind(data.frame(enc = cov_cord %>% codonTable %>% ENC(stop.rm=FALSE),              # Calculate Effective Number of Codons (including STOP codons)
                         cov_cord %>% letterFrequency("GC",as.prob=TRUE)*100 %>% as.vector, # Calculate % GC content
                         cov_cord %>% letterFrequency(c("A","C","G","T")),                  # Nucleotide counts
                         cov_cord %>% dinucleotideFrequency,                                # Dinucleotide counts
                         cov_cord %>% DNAStringSet(start=1) %>% 
                           dinucleotideFrequency(step = 3) %>% as.data.frame %>%
                           rename_all(., ~ paste0(., "_p1")),                               # Dinucleotide counts between positions 1-2 only
                         cov_cord %>% DNAStringSet(start=2) %>% 
                           dinucleotideFrequency(step = 3) %>% as.data.frame %>%
                           rename_all(., ~ paste0(., "_p2")),                               # Dinucleotide counts between positions 2-3 only
                         cov_cord %>% DNAStringSet(start=3) %>% 
                           dinucleotideFrequency(step = 3) %>% as.data.frame %>%
                           rename_all(., ~ paste0(., "_p3")),                               # Dinucleotide counts between positions 3-1 only
                         cov_cord %>% codonTable %>% codonCounts                            # Codon counts
                         #                        cov_cord %>% oligonucleotideFrequency(6, step=3)                   # Codon pair counts - NOT USING FOR NOW
)) %>% rename_at(vars(G.C), ~ "GC_content")

## Helper functions to check correct adjustment of CDS start locations
# dodgyaccessions <- merge(cov_wg_df, allcov_meta_df, all.x = TRUE, by = "accessionversion") %>% 
#   mutate(lengthdiff = CDS_length - meta_length) %>% filter(abs(lengthdiff) > 1000) %>% arrange(-abs(lengthdiff))
# merge(cov_wg_df, allcov_meta_df, all.x = TRUE, by = "accessionversion") %>%
#   filter(length > 10000) %>%
#   ggplot(aes(x = meta_length, y = CDS_length)) +
#   # scale_x_continuous(limits = c(10000,50000)) +
#   # scale_y_continuous(limits = c(10000, 50000)) +
#   geom_point()

cov_enc_df %>% filter(is.na(protein)) %>% nrow # n missing protein annotation
cov_enc_df %>% filter(is.na(gene)) %>% nrow # n missing gene annotation

# Clean up metadata and assign each cds as whole S, S1, S2 or other, and highlight problematic cds
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", protein, ignore.case=TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", protein, ignore.case=TRUE) ~ "S2",
  grepl("spike|^surface|surface gly|s gly|s prot|peplom|protein S$|^S$", protein, ignore.case=TRUE)  ~ "S")
) %>% replace_na(list(seqtype = "other"))

# Override with gene information where available
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", gene, ignore.case=TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", gene, ignore.case=TRUE) ~ "S2",
  TRUE ~ seqtype))

# Flag truncated or mutant proteins to exclude
cov_enc_df %<>% mutate(include = case_when(grepl("trunc|mutant", protein, ignore.case=TRUE) ~ "N")) %>%
  replace_na(list(include = "Y"))

## Code to assist choosing regexp for assigning seqtype
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 0) %>% kable
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 100) %>% arrange(Freq) %>% kable #most common protein labels not captured
# allcov_meta_df %>% filter(as.character(accessionversion) %in% (cov_enc_df %>% filter(is.na(protein)) %>% pull(accession) %>% unique)) %>% nrow # number of covs with unlabelled spike proteins

# Add in taxids corresponding to accesion nos - plyr::join is used to retain row order
cov_enc_df <- plyr::join(allcov_meta_df %>% select(accessionversion, taxid, complete),
                         cov_enc_df,
                         by = "accessionversion",
                         type = "right")

# Assign partial if labelled as complete spike but too small to be complete spike
cov_enc_df %<>% mutate(complete = ifelse(seqtype == "S" & length < min_spike_length, "partial_spike", complete))

# Calculate count, mean length/ENC/GC content of cds for S, S1, S2, other per coronavirus taxid, including only complete cds and merge with main dataset
allcov_df <- merge(allcov_df,
                   cbind(cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           dplyr::count(seqtype) %>% 
                           spread(seqtype, n, fill=0) %>% 
                           rename_at(vars(-taxid), ~ paste0("n_",.)) %>%
                           as.data.frame,
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_length = mean(length, na.rm=TRUE)) %>% 
                           spread(seqtype, mean_length) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_length_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)],
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_enc = mean(enc, na.rm=TRUE)) %>% 
                           spread(seqtype, mean_enc) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_enc_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)],
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_GC = mean(GC_content, na.rm=TRUE)) %>% 
                           spread(seqtype, mean_GC) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_GC_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)]),
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

# Calculate nuc, dinuc, codon counts, GC content of whole genome by summing over all proteins in each sequence
cov_wg_df <- cov_enc_df %>% 
  filter(include == "Y" & complete == "whole_genome") %>% 
  select(-include, -seqtype, -length, -protein, -gene, -title, -complete, -loc_1, -loc_2, -loc_3, -loc_4, -start_new) %>%
  group_by(accessionversion) %>% 
  mutate_if(is.numeric, sum) %>%
  distinct() %>%
  ungroup %>%
  mutate(GC_content = 100*(G+C)/(A+C+G+T), CDS_length = A+G+C+T) %>%
  filter(CDS_length > min_wgs_length, CDS_length < max_wgs_length)   # Impose length restrictions to filter out mislablled data

# Calculate ENC of whole genome
cov_wg_df$enc <- cov_wg_df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% codonTable() %>% ENC

# Calculate mean CDS length, ENC, GC content of whole genomes per coronavirus taxid and merge with main dataset
allcov_df <- merge(allcov_df,
                   cov_wg_df %>%
                     group_by(taxid) %>% 
                     summarise(mean_wgs_length = mean(CDS_length, na.rm=TRUE),
                               mean_wgs_enc = mean(enc, na.rm=TRUE),
                               mean_wgs_GC = mean(GC_content, na.rm=TRUE)),
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

# Is it better to plot length/GC/ENC of ALL whole genomes than mean within taxid?



# ### why is mean_wgs_enc != mean_wgs_length in terms of NAs? it's those missing sequences again - probably want to investigate this a little as it is 78 viruses!
# hi <- allcov_meta_df %>% filter(complete == "whole_genome") %>% distinct(accessionversion) %>% pull(accessionversion)
# # 2954
# lo <- cov_wg_df %>% filter(accessionversion %in% (allcov_meta_df %>% filter(complete == "whole_genome") %>% pull(accessionversion))) %>% pull(accessionversion)
# # 2876
# missingaccession <- hi[!(hi %in% lo)]
# missinguid <- allcov_meta_df %>% filter(accessionversion %in% missingaccession) %>% pull(uid)
# missinguid %>% entrez_fetch(db = "nuccore",id = ., rettype="fasta_cds_na")
# save(missingaccession, missinguid, file = "data//missingaccession.RData")
# 
# # length is pulled from metadata,  whereas GC/ENC must be taken from sequence, so even if it returns "/n" we can still get length
# allcov_df %>% filter(!is.na(mean_wgs_length)) %>% nrow=
# allcov_df %>% filter(!is.na(mean_wgs_GC)) %>% nrow
# 
# missingids <- allcov_df %>% filter(!is.na(mean_wgs_length) & is.na(mean_wgs_GC)) %>% pull(childtaxa_id)
# missingaccession <- allcov_meta_df %>% filter(taxid %in% missingids & complete == "whole_genome") %>% pull(accessionversion)
# allcov_meta_df %>% filter(accessionversion %in% missingaccession) %>% pull(uid)

# just pulling in "\n"...  THE REASON IS THE SEQUENCES ARE NOT ACTUALLY ANNOTATED WITH CODING SEQUENCE REGIONS AND NONCODING SEQUENCE REGIONS
# IF WE REALLY NEED THESE I SUPPOSE WE COULD TRY ALIGNING?
# CAN THEY BE FILTERED EARLIER TO KEEP NUMBERS UPDATED, E.G. IN ENTREZ_SEARCH AND ENTREZ_SUMMARY

#################################
# Calculate values of interest  #
#################################

# Count total sequences for various criteria
n_cov <- nrow(allcov_df) # n coronaviruses
n_cov_seq <- allcov_df %>% filter(n_seqs > 0) %>% nrow # n coronaviruses with at least 1 sequence
n_cov_seq_eid2 <- allcov_df %>% filter(n_seqs > 0 & eid2_n_species > 0) %>% nrow # n coronaviruses with at least 1 sequence and at least 1 EID2 host
n_seq <- allcov_meta_df %>% nrow # n nucleotide sequence entries
n_cds <- cov_enc_df %>% nrow # n coding sequences

# # 492 viruses should have complete spike protein based on metadata, but only 420 do based on the cds labels in the FASTA!miss_spike <- tab492[!(tab492 %in% tab420)]
# tab420 <- allcov_df %>% filter(!is.na(n_S) & n_S > 0) %>% pull(childtaxa_id)
# tab492 <-  allcov_meta_df %>% group_by(taxid, complete) %>% count(complete) %>% spread(complete, n, fill=0) %>% rename_at(vars(-taxid), ~ paste0("n_",.)) %>% filter(n_complete_spike > 0 | n_whole_genome > 0) %>% pull(taxid) %>% as.character
# 
# miss_spike <- tab492[!(tab492 %in% tab420)] # Identify taxids missing the spike protein
# 
# allcov_df %>% filter(childtaxa_id %in% miss_spike) %>% pull(n_other) %>% table(exclude = FALSE) # 8 come through in other sequences, but 64 have NAs
# =
# # For the 8, metadata describes them as complete genome...
# 
# allcov_meta_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & !(is.na(n_other))) %>% pull(childtaxa_id))) %>% pull(title)
# 
# # ...yet they don't have any spike proteins (or none annotated)!
# 
# cov_enc_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & !(is.na(n_other))) %>% pull(childtaxa_id))) %>% pull(protein)
# 
# # For the 64 - entrez_fetch just.. fails to fetch them and returns blanks???
# 
# allcov_meta_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & is.na(n_other)) %>% pull(childtaxa_id))) %>% pull(uid) %>% as.character %>% entrez_fetch(db = "nuccore",id = ., rettype="fasta_cds_na")
# 
# # So 420 is the correct number of viruses with a spike protein...
# sum(allcov_df$n_S, na.rm=TRUE) # 5213 individual spike protein sequences

###################################
# Establish spike protein dataset #
###################################

# Filter dataset to final complete spike protein dataset
cov_spikes_df <- cov_enc_df %>% filter(seqtype == "S" & include == "Y" & complete != "partial_spike")

# Write FASTA of final complete spike protein dataset
cov_cord[which(cov_enc_df$seqtype == "S" & cov_enc_df$include == "Y" & cov_enc_df$complete != "partial_spike"), ] %>%
  writeXStringSet(., filepath = "data\\allcov_spikes.fasta", format="fasta")

write(unlist(allcov_FASTA), file="data\\allcov_GenBank.fasta")

merge(merge(cov_spikes_df %>% select(accessionversion, taxid, title, complete, length),
            allcov_meta_df %>% select(accessionversion, createdate, subname),
            by = "accessionversion",
            all.x = TRUE),
      allcov_df %>% select(childtaxa_id, childtaxa_name, genus),
      by.x = "taxid",
      by.y = "childtaxa_id",
      all.x = TRUE) %>% select(accessionversion, title, taxid, genus, childtaxa_name, createdate, complete, length, subname) %>% write.csv(., "data\\allcov_spikes_metadata.csv")


#######################################
# Calculate genome composition biases #
#######################################

cov_spikes_df %<>% calc_composition_bias
cov_wg_df %<>% calc_composition_bias

# Append genus and virus name information
cov_spikes_df %<>% left_join(allcov_df %>% select(childtaxa_id, childtaxa_name, genus) %>% mutate(genus = gsub("0","unclassified",genus)),
                             by = c("taxid" = "childtaxa_id"))

cov_wg_df %<>% left_join(allcov_df %>% select(childtaxa_id, childtaxa_name, genus) %>% mutate(genus = gsub("0","unclassified",genus)),
                         by = c("taxid" = "childtaxa_id"))

# Append host metadata information
cov_spikes_df %<>% left_join(allcov_meta_df %>% select(accessionversion, class_name, order_name, family_name, genus_name, species_name),
                             by = "accessionversion")

cov_wg_df %<>% left_join(allcov_meta_df %>% select(accessionversion, class_name, order_name, family_name, genus_name, species_name),
                         by = "accessionversion")





############
# Heatmaps #
############

# Take a random sample just to illustrate heatmap size for 420 viruses
mat_data <- cov_spikes_df %>% sample_n(420) 

heatmap.2(mat_data %>% 
            set_rownames(.$accessionversion) %>% 
            select(matches("^[A|C|G|T][A|C|G|T]_Bias$")) %>% 
            as.matrix,
          density.info="none", trace="none", margins=c(12,12), dendrogram="row", Colv="NA", 
          # labRow = cov_pca_df$taxid,   # specify row labels = species?
          RowSideColors = mat_data %>% mutate(rowsidecol = case_when(                       # Set side colours using same genus colours as ggplots elsewhere
            genus == "Alphacoronavirus" ~ "#F8766D",
            genus == "Betacoronavirus" ~ "#A3A500",
            genus == "Deltacoronavirus" ~ "#00BF7D",
            genus == "Gammacoronavirus" ~ "#00B0F6",
            genus == "unclassified" ~ "#E76BF3",
          )) %>% pull(rowsidecol),
          col = colorRampPalette(c("dodgerblue", "gray10", "firebrick2"))(n = 55),
          breaks = c(seq(0,0.95,length=25), seq(0.951,1.05,length=6), seq(1.051,2.0,length=25)))

par(lend = 1)           # square line ends for the color legend
legend(x=0.75,y=0.2, legend = mat_data$genus %>% unique %>% sort, fill = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"), ncol = 5)










#################################################
# OLD TEMP place to check slippage/frameshifts  #
#################################################

# ### MANUAL CORRECTION HELPER CODE
# seq <- f$`lcl|NC_039196.1_cds_YP_009512959.1_3` %>% as.vector() %>% as.list %>% do.call(paste0, .)
# substring(seq,seq(1,nchar(seq),3),seq(3,nchar(seq),3))
# aminoseq <- "MSSHQIQQVKHGLESLQEIKNNPPSSQDVNLAREIYESIRQTGTSSVQGGAIAGDNITSGGNNDSMYSQGPSPPISSVNKNIEGPTGFDHSGLWDPEGNLCMLFESDDDENHYSEINGRSSAIEGLDEQDNENSIIKQPGNQCTEGVSKTDSSLSSQETTLSVGGSDIPGAGISTCASLDITVNELEDATVRNSNNMKGNWPIPKLLVKPPPRVKTSVDHSNPLKGGHRREISLTWDGDYIIREEWCNPICTPIYSTCKRLQCRCKQCPSTCPKCE"
# str_locate_all(aminoseq, "K")
# substring(seq,seq(1,nchar(seq),3),seq(3,nchar(seq),3))[225:227]
# substring(seq, (138)*3, (139)*3-1)
# 
# n <- 225
# substring(seq, (n-1)*3+1, (n)*3) # Gives nth codon
# # If n is the codon the insert goes after, insert starts from position n*3+1
# 
# ins_pos <- 421 # Insertion position
# ins_seq <- "G" # Insertion sequence
# gsub(paste0('^([a-z]{',ins_pos-1,'})([a-z]+)$'), paste0('\\1',ins_seq,'\\2'), seq, perl=TRUE) %>% writeClipboard
