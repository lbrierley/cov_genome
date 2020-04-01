###############################
# Extract sequences from NCBI #
###############################

# Extract all viruses in family Coronaviridae and their taxonomic IDs
allcov_df <- downstream("Coronaviridae", db="ncbi")$Coronaviridae %>% filter(rank %in% c("no rank", "species"))

# Find and obtain results based on taxonomy IDs of all coronaviruses
if (load_prev_seqs == TRUE){
  load(file="allcov_results.RData")
} else {
  allcov_results <- lapply(as.list(allcov_df$childtaxa_id), function(x) Seq_searcher(x))
  save(allcov_results, file="allcov_results.RData")
}

# Append sequence counts of viruses found
allcov_df <- data.frame(allcov_df, n_seqs = lapply(allcov_results, function(x) x$count) %>% unlist())

# Extract IDs of individual sequences
allcov_seqids <- lapply(allcov_results, function(x) x$ids) %>% Filter(length, .) 

# Fetch metadata information for those sequences (title, accession). NB entrez_summary doesn't accept empty objects, so filtering out those IDs with no sequences
if (load_prev_seqs == TRUE){
  load(file="allcov_found_summaries.RData")
} else {
  allcov_found_summaries <- pblapply(allcov_seqids, function(x) Seq_summary(x))
  save(allcov_found_summaries, file="allcov_found_summaries.RData")
}

# Fetch sequences as FASTA
if (load_prev_seqs == TRUE){
  load(file="allcov_FASTA.RData")
} else {
  allcov_FASTA <- pblapply(allcov_seqids, function(x) Seq_FASTA(x))
  save(allcov_FASTA, file="allcov_FASTA.RData")
}

# Write FASTA file
write(unlist(allcov_FASTA), file="allcov_GenBank.fasta")

# Identify n sequences for SARS-CoV-2
allcov_df %>% filter(childtaxa_id == 2697049)

# NB numbers don't match between this and SARS-CoV-2 special page on GenBank because some of them have incorrectly just been filed under taxid 694009 for SARS-CoV parent species...
# Load in FASTA from manually downloaded GISAID repository (2020-03-27 18:15:04 GMT)
sarscov2_GenBank_FASTA <- read.fasta("sarscov2_GenBank.fasta")
sarscov2_GISAID_FASTA <- read.fasta("sarscov2_GISAID.fasta")

length(sarscov2_GenBank_FASTA) # total CDS
sarscov2_GenBank_FASTA[lapply(sarscov2_GenBank_FASTA, function(x) x %>% 
                                attributes %>% 
                                .$Annot %>% 
                                grepl("gene=S", .)) %>% unlist] %>%  length # total S sequences - should be larger but it's smaller?? some might be missing annotation?

length(sarscov2_GISAID_FASTA) # total complete genomes

########################
# Add host information #
########################

eid2_cov <- read.csv("species_species_dataset_aves_mammals.csv", fileEncoding="UTF-8-BOM")  %>% filter(PathogenTaxID %in% allcov_df$childtaxa_id)

# Recast data to create 0/1 outcome variables as to whether each virus is known to infect each host type
eid2_outcomes <- dcast(eid2_cov, PathogenTaxID ~Hc1, length)
eid2_outcomes %<>% mutate_at(vars(-PathogenTaxID), ~1 * (. > 0))
eid2_outcomes %<>% mutate(n_groups = rowSums(eid2_outcomes[,-1]))
eid2_outcomes %<>% mutate_at(vars(-PathogenTaxID, -n_groups), funs(factor))

# Count number of host species infected
# Including the following as host species based on EID2 data structure: humans, domestic subspecies (e.g canis lupus separate to canis lupus familiaris)
eid2_outcomes <- merge(eid2_outcomes,
                       eid2_cov %>% 
                         group_by(PathogenTaxID) %>% 
                         summarise(Hostspp = n_distinct(HostTaxID)),
                       by = "PathogenTaxID")

allcov_df <- merge(allcov_df, eid2_outcomes,
                   by.x = "childtaxa_id",
                   by.y = "PathogenTaxID",
                   all.x = TRUE) %>% replace(., is.na(.), 0)