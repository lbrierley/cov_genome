###############################
# Extract sequences from NCBI #
###############################

# Extract all viruses in family Coronaviridae and their taxonomic IDs
allcov_df <- downstream("Coronaviridae", db = "ncbi")$Coronaviridae %>% filter(rank %in% c("no rank", "species"))

# Find and obtain results based on taxonomy IDs of all coronaviruses
if (load_prev_seqs == TRUE) {
  load(file = "data\\allcov_results.RData")
} else {
  allcov_results <- lapply(as.list(allcov_df$childtaxa_id), function(x) Seq_searcher(x))
  save(allcov_results, file = "data\\allcov_results.RData")
}

# Load IDs of sequences to exclude as they're known a priori to have no CDS annotation
load("data\\missingaccession.RData")

# Append genus name, sequence counts of viruses found
allcov_df <- data.frame(allcov_df,
                        n_seqs = lapply(allcov_results, function(x) {
                          x$ids %>%
                            subset(!(. %in% (missinguid %>% as.character()))) %>%
                            length()
                        }) %>% unlist()
)

genus_list <- lapply(allcov_df$childtaxa_id %>% classification(), function(x) {
  x %>%
    filter(rank == "genus") %>%
    pull(name)
})
is.na(genus_list) <- lengths(genus_list) == 0

allcov_df$genus <- genus_list %>%
  unlist() %>%
  unname()

# Extract IDs of individual sequences
allcov_seqids <- lapply(allcov_results, function(x) x$ids %>% 
                          subset(!(. %in% (missinguid %>% as.character())))) %>% 
  Filter(length, .) # entrez_summary doesn't accept empty objects, so filtering out those taxonomic IDs with no sequences

# Fetch metadata information for those sequences (title, accession). 
if (load_prev_seqs == TRUE) {
  load(file = "data\\allcov_found_summaries.RData")
} else {
  allcov_found_summaries <- pblapply(allcov_seqids, function(x) Seq_summary(x))
  save(allcov_found_summaries, file = "data\\allcov_found_summaries.RData")
}

# Fetch sequences as FASTA
if (load_prev_seqs == TRUE) {
  load(file = "data\\allcov_FASTA.RData")
} else {
  allcov_FASTA <- pblapply(allcov_seqids, function(x) Seq_FASTA(x))
  save(allcov_FASTA, file = "data\\allcov_FASTA.RData")
}

# Write FASTA file
write(unlist(allcov_FASTA), file = "data\\allcov.fasta")

##############################
# Clean and process metadata #
##############################

# Extract metadata into tidy df format
allcov_found_summaries %<>% flatten

class(allcov_found_summaries) <- c("list", "esummary_list") # coerce to esummary list

allcov_meta_df <- sapply(
  c("uid", "caption", "title", "extra", "gi", "createdate", "updatedate", "flags", "taxid", "slen", "biomol", "moltype", "topology", "sourcedb", "segsetsize", "projectid", "genome", "subtype", "subname", "assemblygi", "assemblyacc", "tech", "completeness", "geneticcode", "strand", "organism", "strain", "biosample", "accessionversion"),
  function(x) extract_from_esummary(allcov_found_summaries, x) %>% as.vector()
) %>% as.data.frame()

allcov_meta_df$meta_length <- extract_from_esummary(allcov_found_summaries, "statistics") %>%
  t() %>%
  .[, 2] %>%
  lapply(., function(x) x[1]) %>%
  unlist()

# Clean up metadata and assign each entry as complete or partial genome
allcov_meta_df %<>% mutate(complete = case_when(
  grepl("complete genome", title, ignore.case = TRUE) ~ "whole_genome", # Assign as whole genome if labelled so
  grepl("partial and complete cds", title, ignore.case = TRUE) ~ "partial_spike", # Else assign partial if mixed
  metadata_title_cleaner(title) == "partial" ~ "partial_spike", # Else apply the metadata cleaner
  metadata_title_cleaner(title) == "complete" ~ "complete_spike",
  grepl("partial", title, ignore.case = TRUE) ~ "partial_spike", # If the complex metadata cleaner can't assign it then just do a simple search match
  grepl("complete", title, ignore.case = TRUE) ~ "complete_spike"
)) %>%
  replace_na(list(complete = "partial_spike")) # Else no mention of completeness, then just assume not complete

# Re-format sub-metadata and add in as columns
sub_meta_df <- apply(allcov_meta_df, 1, function(x) {
  x["subname"] %>%
    str_split("\\|") %>%
    unlist() %>%
    t() %>%
    data.frame() %>%
    set_colnames(x["subtype"] %>%
                   str_split("\\|") %>%
                   unlist())
}) %>% bind_rows()

sub_meta_df <- sub_meta_df[, which(colnames(sub_meta_df) != "")] # Remove blank columns created by empty metadata
sub_meta_df %<>% rename_all(function(x) {
  paste0("meta_", x)
}) # Rename to label them as sub metadata
allcov_meta_df %<>% cbind(sub_meta_df)
rm(sub_meta_df)

# Wipe host metadata for non-natural or recombinant viruses missed by previous filters in order to exclude these
allcov_meta_df %<>%
  mutate(
    meta_host =
      case_when(
        grepl("recomb|vacc", subname) ~ NA_character_,
        TRUE ~ meta_host
      )
  )

if (load_prev_seqs == TRUE) {
  load(file = "data\\tax_matcher.RData") # Load previous host name cross-referencing
} else {
  
  # Setup empty list
  list_tax_matcher <- vector("list", nrow(allcov_meta_df))
  
  # Resolve host metadata to taxonomic species
  for (i in 1:nrow(allcov_meta_df)) {
    list_tax_matcher[[i]] <-
      if (!(allcov_meta_df %>% slice(i) %>% pull(meta_host) %>% is.na()) == T) {
        allcov_meta_df %>%
          slice(i) %>%
          pull(meta_host) %>%
          gsub(" \\(.*$|\\;.*$|\\,.*$", "", .) %>%
          gsub(" sp\\.| spp\\.", "", .) %>%
          name2taxid(out_type = "summary") %>% # If host metadata is present, resolve and add in accessionversion
          mutate(accessionversion = allcov_meta_df %>%
                   slice(i) %>%
                   pull(accessionversion))
      } else {
        tibble(
          name_txt = character(),
          tax_id = character(),
          accessionversion = factor()
        ) # If host metadata is not present, return an empty tibble so list can still be collapsed with rbindlist
      }
    
    if (list_tax_matcher[[i]] %>% nrow() == 0) { # If empty (either no host metadata, or host metadata not found), return given term alongside NA, and accessionversion
      list_tax_matcher[[i]] <- tibble(
        name_txt = allcov_meta_df %>%
          slice(i) %>%
          pull(meta_host),
        tax_id = NA,
        accessionversion = allcov_meta_df %>%
          slice(i) %>%
          pull(accessionversion)
      )
    }
  }
  
  # Collapse to data frame
  df_tax_matcher <- data.table::rbindlist(list_tax_matcher)
  
  # Manually resolve remaining unmatched host labels to lowest possible unambiguous taxonomic level
  df_tax_matcher %<>% mutate(
    tax_id = case_when(
      
      # Fix common names that match multiple taxonomic ids with their lowest common parent taxid
      name_txt == "Camel" | name_txt == "camel" ~ "9836", 
      name_txt == "mink" ~ "169418",
      name_txt == "rat" ~ "10114",
      name_txt == "mouse" ~ "10088",
      !is.na(tax_id) ~ tax_id,
      
      # Fix common domestic species
      grepl("Gallus gallus|poultry|chicken|chikcen|breeder|layer|broiler|\\<hen\\>|\\<hens\\>|chick with 18 days|embryonated egg", name_txt, ignore.case = TRUE) ~ "9031", 
      grepl("Bos taurus|cow|calf", name_txt, ignore.case = TRUE) ~ "9913",
      grepl("Sus scrofa|pig|porcine", name_txt, ignore.case = TRUE) ~ "9823",
      grepl("Equus caballus|horse", name_txt, ignore.case = TRUE) ~ "9796",
      grepl("Felis catus|feline", name_txt, ignore.case = TRUE) ~ "9685",
      grepl("Canis lupus familiaris|canine|dog", name_txt, ignore.case = TRUE) ~ "9615",
      grepl("Equus africanus asinus|donkey", name_txt, ignore.case = TRUE) ~ "9793",
      grepl("Rhinilophus ferrumequinum", name_txt, ignore.case = TRUE) ~ "59479", 
      
      # Fix clear instances of errors or undetectable otherwise unambiguous species
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
      
      #Else keep the automated taxonomic id
      TRUE ~ tax_id
    )
  ) %>%
    distinct()
  
  # For each host taxonomic ID retrieve and arrange into data frame class, order, family, genus, species
  host_taxize_full_class <- df_tax_matcher %>%
    filter(!is.na(tax_id)) %>%
    pull(tax_id) %>%
    unique() %>%
    classification(db = "ncbi") %>%
    lapply(function(x) {
      x %>%
        data.frame() %>%
        filter(rank == "class" | rank == "order" | rank == "family" | rank == "genus" | rank == "species") %>%
        # Reshape taxonomy tables to give one row per host with family/genus/species*id/name columns
        melt(id = "rank") %>%
        unite(rank_var, c("rank", "variable")) %>%
        spread(rank_var, value)
    }) %>%
    bind_rows() %>%
    mutate(tax_id = df_tax_matcher %>% filter(!is.na(tax_id)) %>% pull(tax_id) %>% unique())
  
  # Merge in full taxonomic information
  df_tax_matcher %<>% left_join(host_taxize_full_class, by = "tax_id")
  
  save(df_tax_matcher, file = "data\\tax_matcher.RData")
}

# Append full metadata host label taxonomy to sequence-level data
allcov_meta_df %<>% left_join(df_tax_matcher, by = "accessionversion")

# Create functional host group variable for use as an ML outcome
allcov_meta_df %<>%
  mutate(
    group_name =
      case_when(
        species_name == "Homo sapiens" ~ "human",
        family_name == "Suidae" ~ "swine",
        family_name == "Camelidae" ~ "camel",
        family_name %in% c("Craseonycteridae", "Hipposideridae", "Megadermatidae", "Pteropodidae", "Rhinolophidae", "Rhinopomatidae") ~ "yinbat",
        family_name %in% c("Emballonuridae", "Furipteridae", "Miniopteridae", "Molossidae", "Mormoopidae", "Mystacinidae", "Myzopodidae", "Natalidae", "Noctilionidae", "Nycteridae", "Phyllostomidae", "Thyropteridae", "Cistugidae", "Vespertilionidae")
        ~ "yangbat",
        order_name == "Rodentia" ~ "rodent",
        order_name == "Carnivora" ~ "carnivore",
        class_name == "Aves" ~ "aves"
      )
  )

# Count host categories per coronavirus and append to species-level data
allcov_df %<>% 
  left_join(allcov_meta_df %>%
              group_by(taxid) %>%
              summarise_at(c("class_name", "order_name", "family_name", "genus_name", "species_name", "group_name"), n_distinct, na.rm = TRUE) %>%
              rename_at(vars(-taxid), ~ paste0("n_", gsub("_name", "", .))),
            by = c("childtaxa_id" = "taxid")
  )

# Read FASTA back in in coRdon format and calculate codon table
cov_cord <- readSet(file = "data\\allcov.fasta")

# Construct summary df for individual CDS
cov_enc_df <- data.frame(
  title = cov_cord %>% names(),
  accessionversion = cov_cord %>% names() %>% str_match("lcl\\|(.*?)_cds_") %>% .[, 2] %>% as.character(),
  gene = cov_cord %>% names() %>% str_match("\\[gene=(.*?)\\]") %>% .[, 2],
  protein = cov_cord %>% names() %>% str_match("\\[protein=(.*?)\\]") %>% .[, 2],
  length = cov_cord %>% width()
)

# Extract and append CDS location information by text extraction from title
suppressWarnings(
  cov_enc_df %<>% cbind(cov_enc_df %>%
                          pull(title) %>%
                          str_match("\\[location=(.*?)\\]") %>%
                          .[, 2] %>%
                          gsub("join\\(|\\)|>|<", "", .) %>%
                          gsub(",", "..", .) %>%
                          str_split(., "\\..") %>%
                          lapply(., function(x) {
                            x %>%
                              t() %>%
                              as.data.frame()
                          }) %>%
                          bind_rows() %>%
                          mutate_all(., as.numeric) %>%
                          rename_all(., ~ gsub("V", "loc_", .)))
)

if(frameshift_correct == TRUE) {
  
  # Correct frameshift/slippage causing double accounting of elements of some ORF1ab CDS by specifying new location to start calculating metrics at
  cov_enc_df %<>% left_join(
    cov_enc_df %>% filter( # Only apply correction where needed
      accessionversion %in% (cov_enc_df %>% filter(grepl("join", title)) %>% pull(accessionversion)) & # Filter to sequences containing a CDS that join overlapping components
        accessionversion %in% (allcov_meta_df %>% filter(complete == "whole_genome") %>% pull(accessionversion)) # Filter to sequences that describe a whole genome
    ) %>% group_by(accessionversion, loc_1) %>% arrange(accessionversion, loc_1) %>% filter(n() > 1) %>% # Filter to sequences with two different CDS that both start at the same location
      mutate(start_new = ifelse(is.na(loc_3), lag(loc_2) - (loc_1 - 1) + 1, 1)) %>% # For the shorter sequence of ORF1a, set new start location based on slippage location of ORF1b
      ungroup() %>%
      select(title, start_new),
    by = "title"
  ) %>%
    replace_na(list(start_new = 1))
  
  # Apply new CDS start locations
  cov_cord %<>% subseq(start = cov_enc_df$start_new)
  
}

# Calculate genomic composition counts for individual CDS
cov_enc_df %<>% cbind(data.frame(
  enc = cov_cord %>% codonTable() %>% ENC(stop.rm = FALSE), # Calculate Effective Number of Codons (including STOP codons)
  cov_cord %>% letterFrequency("GC", as.prob = TRUE) * 100 %>% as.vector(), # Calculate % GC content
  cov_cord %>% letterFrequency(c("A", "C", "G", "T")), # Nucleotide counts
  cov_cord %>% dinucleotideFrequency(), # Dinucleotide counts
  cov_cord %>%
    DNAStringSet(start = 1) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p1")), # Dinucleotide counts between positions 1-2 only
  cov_cord %>%
    DNAStringSet(start = 2) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p2")), # Dinucleotide counts between positions 2-3 only
  cov_cord %>%
    DNAStringSet(start = 3) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p3")), # Dinucleotide counts between positions 3-1 only
  cov_cord %>%
    codonTable() %>%
    codonCounts() # Codon counts
)) %>% rename_at(vars(G.C), ~"GC_content")

cov_enc_df %>%
  filter(is.na(protein)) %>%
  nrow() # n missing protein annotation
cov_enc_df %>%
  filter(is.na(gene)) %>%
  nrow() # n missing gene annotation

# Clean up metadata and assign each cds as whole S, S1, S2 or other
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", protein, ignore.case = TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", protein, ignore.case = TRUE) ~ "S2",
  grepl("spike|^surface|surface gly|s gly|s prot|peplom|protein S$|^S$", protein, ignore.case = TRUE) ~ "S"
)) %>% replace_na(list(seqtype = "other"))

# Override with gene information where available
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", gene, ignore.case = TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", gene, ignore.case = TRUE) ~ "S2",
  TRUE ~ seqtype
))

# Flag truncated or mutant proteins to exclude
cov_enc_df %<>% mutate(include = case_when(grepl("trunc|mutant", protein, ignore.case = TRUE) ~ "N")) %>%
  replace_na(list(include = "Y"))

# Add in taxids corresponding to accesion nos - plyr::join is used to retain row order
cov_enc_df <- plyr::join(allcov_meta_df %>% select(accessionversion, taxid, complete),
                         cov_enc_df,
                         by = "accessionversion",
                         type = "right"
)

# Assign partial if labelled as complete spike but too small to be complete spike
cov_enc_df %<>% mutate(complete = case_when(
  seqtype == "S" & length < min_spike_length ~ "partial_spike",
  TRUE ~ complete))

cov_wg_df <- cov_enc_df %>%
  select(-seqtype, -length, -protein, -gene, -title, -loc_1, -loc_2, -loc_3, -loc_4)

if(frameshift_correct == TRUE) {
  cov_wg_df %<>% select(-start_new)
}

# Calculate nuc, dinuc, codon counts, GC content of whole genome by summing over all proteins in each sequence
cov_wg_df %<>%
  filter(include == "Y" & complete == "whole_genome") %>%
  select(-include, -complete) %>%
  group_by(accessionversion) %>%
  mutate_if(is.numeric, sum) %>%
  distinct() %>%
  ungroup() %>%
  mutate(GC_content = 100 * (G + C) / (A + C + G + T), CDS_length = A + G + C + T) %>%
  filter(CDS_length > min_wgs_length) # Impose length restrictions to filter out mislabelled data

if(frameshift_correct == TRUE) {
  cov_wg_df %<>% filter(CDS_length < max_wgs_length) # Only impose upper length restriction if correcting for frameshift (as if not, sizes as-read will be >> size as-is)
}

# Calculate ENC of whole genome
cov_wg_df$enc <- cov_wg_df %>%
  select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>%
  codonTable() %>%
  ENC()

# Filter dataset to final complete spike protein dataset
cov_spikes_df <- cov_enc_df %>% filter(seqtype == "S" & include == "Y" & complete != "partial_spike")

#######################################
# Calculate genome composition biases #
#######################################

cov_spikes_df %<>% calc_composition_bias
cov_wg_df %<>% calc_composition_bias

# Append genus and virus name information
cov_spikes_df %<>% left_join(allcov_df %>% select(childtaxa_id, childtaxa_name, genus) %>% mutate(genus = gsub("0", "unclassified", genus)),
                             by = c("taxid" = "childtaxa_id")
)

cov_wg_df %<>% left_join(allcov_df %>% select(childtaxa_id, childtaxa_name, genus) %>% mutate(genus = gsub("0", "unclassified", genus)),
                         by = c("taxid" = "childtaxa_id")
)

# Append host metadata information
cov_spikes_df %<>% left_join(allcov_meta_df %>% select(accessionversion, class_name, order_name, family_name, genus_name, species_name, group_name),
                             by = "accessionversion"
)

cov_wg_df %<>% left_join(allcov_meta_df %>% select(accessionversion, class_name, order_name, family_name, genus_name, species_name, group_name),
                         by = "accessionversion"
)

#############
# Save data #
#############

save(allcov_df,
     cov_spikes_df,
     cov_wg_df,
     file = paste0("data\\cov_ML_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData")
)