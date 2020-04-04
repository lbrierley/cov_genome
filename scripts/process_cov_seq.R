###############################
# Extract sequences from NCBI #
###############################

# ALL CORONAVIRIDAE

# Extract all viruses in family Coronaviridae and their taxonomic IDs
allcov_df <- downstream("Coronaviridae", db="ncbi")$Coronaviridae %>% filter(rank %in% c("no rank", "species"))

# Find and obtain results based on taxonomy IDs of all coronaviruses
if (load_prev_seqs == TRUE){
  load(file="data\\allcov_results.RData")
} else {
  allcov_results <- lapply(as.list(allcov_df$childtaxa_id), function(x) Seq_searcher(x))
  save(allcov_results, file="data\\allcov_results.RData")
}

# Append sequence counts of viruses found
allcov_df <- data.frame(allcov_df, n_seqs = lapply(allcov_results, function(x) x$count) %>% unlist())

# Extract IDs of individual sequences
allcov_seqids <- lapply(allcov_results, function(x) x$ids) %>% Filter(length, .)

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

##############################
# Clean and process metadata #
##############################

# Extract metadata into df format
allcov_found_summaries %<>% flatten
class(allcov_found_summaries) <-  c("list", "esummary_list") # coerce to esummary list
allcov_meta_df <- sapply(c("uid", "caption", "title", "extra", "gi", "createdate", "updatedate", "flags", "taxid", "slen", "biomol", "moltype", "topology", "sourcedb", "segsetsize", "projectid", "genome", "subtype", "subname", "assemblygi", "assemblyacc", "tech", "completeness", "geneticcode", "strand", "organism", "strain", "biosample", "accessionversion")
                         , function(x) extract_from_esummary(allcov_found_summaries, x) %>% as.vector) %>% as.data.frame
allcov_meta_df$length <- extract_from_esummary(allcov_found_summaries, "statistics") %>% t %>% .[,2] %>% lapply(., function(x) x[1]) %>% unlist

# Clean up metadata and assign each entry as complete or partial genome
allcov_meta_df %<>% mutate(complete = case_when(
  grepl("complete genome", title, ignore.case=TRUE) ~ "whole_genome",    # Assign as whole genome if labelled so
  metadata_title_cleaner(title) == "partial" ~ "partial_spike",          # Else apply the metadata cleaner
  metadata_title_cleaner(title) == "complete"  ~ "complete_spike",
  grepl("partial", title, ignore.case=TRUE) ~ "partial_spike",           # If the complex metadata cleaner can't assign it (which is only 3% of records) then just do a simple search match
  grepl("complete", title, ignore.case=TRUE) ~ "complete_spike")) %>% 
 replace_na(list(complete = "complete_spike"))                           # Else no mention of completeness (which is only 0.3% of records), then just assume complete spike

## Code to assist choosing regexp for assigning complete
# table(allcov_meta_df$completeness, exclude = NULL) #MOST NOT LABELLED
# allcov_meta_df %>% filter(complete == "partial_spike" & completeness == "complete") %>% select(title) # check against metadata values where they exist
# allcov_meta_df %>% filter(complete == "complete_spike" & completeness == "partial") %>% select(title) # looks like the metadata "completeness" field is wrong or at least conflicts with title

# Read FASTA back in in coRdon format and calculate codon table
cov_cord <- readSet(file="data\\allcov_GenBank.fasta")

# Construct summary df for individual cds
cov_enc_df <- data.frame(title = cov_cord %>% names(),
                         accession = cov_cord %>% names() %>% str_match("lcl\\|(.*?)_cds_") %>% .[,2] %>% as.character,
                         gene = cov_cord %>% names() %>% str_match("\\[gene=(.*?)\\]") %>% .[,2],
                         protein = cov_cord %>% names() %>% str_match("\\[protein=(.*?)\\]") %>% .[,2],
                         length = cov_cord %>% width,
                         enc = cov_cord %>% codonTable %>% ENC(stop.rm=FALSE)) %>% # Calculate Effective Number of Codons (including STOP codons)
  cbind(alphabetFrequency(cov_cord, as.prob = TRUE)[,1:4])

cov_enc_df %>% filter(is.na(protein)) %>% nrow # n missing protein annotation
cov_enc_df %>% filter(is.na(gene)) %>% nrow # n missing gene annotation

# Clean up metadata and assign each cds as whole S, S1, S2 or other, and highlight problematic cds
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", protein, ignore.case=TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", protein, ignore.case=TRUE) ~ "S2",
  grepl("spike|^surface|surface gly|s gly|s prot|peplom|protein S$|^S$", protein, ignore.case=TRUE)  ~ "S")
) %>% replace_na(list(seqtype = "other"))

cov_enc_df %<>% mutate(include = case_when(grepl("trunc|mutant", protein, ignore.case=TRUE) ~ "N")) %>%
  replace_na(list(include = "Y"))

## Code to assist choosing regexp for assigning seqtype
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 0) %>% kable
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 100) %>% arrange(Freq) %>% kable #most common protein labels not captured
# allcov_meta_df %>% filter(as.character(accessionversion) %in% (cov_enc_df %>% filter(is.na(protein)) %>% .$accession %>% unique)) %>% nrow # number of covs with unlabelled spike proteins

# Add in taxids corresponding to accesion nos
cov_enc_df <- merge(allcov_meta_df %>% select(accessionversion, taxid, complete),
                    cov_enc_df,
                    by.x = "accessionversion",
                    by.y = "accession",
                    all.y = TRUE)

# Calculate count, mean length/ENC/GC content of cds for S, S1, S2, other per coronavirus taxid, including only complete cds and merge with main dataset

allcov_df <- merge(allcov_df,
                   cbind(cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           count(seqtype) %>% 
                           spread(seqtype, n, fill=0) %>% 
                           rename_at(vars(-taxid), ~ paste0("n_",.)) %>%
                           as.data.frame,
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_length = mean(length)) %>% 
                           spread(seqtype, mean_length) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_length_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)],
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_enc = mean(enc)) %>% 
                           spread(seqtype, mean_enc) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_enc_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)],
                         cov_enc_df %>% 
                           filter(include == "Y" & complete != "partial_spike") %>% 
                           group_by(taxid, seqtype) %>% 
                           summarise(mean_GC = mean(G+C)) %>% 
                           spread(seqtype, mean_GC) %>% 
                           rename_at(vars(-taxid), ~ paste0("mean_GC_",.)) %>% 
                           as.data.frame %>%
                           .[c(2:5)]),
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

# Calculate mean length of complete genomes per coronavirus taxid and merge with main dataset

allcov_df <- merge(allcov_df,
                   allcov_meta_df %>% 
                     filter(complete == "whole_genome") %>% 
                     group_by(taxid) %>% 
                     summarise(mean_wgs_length = mean(length)) %>% 
                     as.data.frame,
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

#################################
# Calculate values of interest  #
#################################

# Count total sequences for various criteria
n_cov <- nrow(allcov_df) # n coronaviruses
n_cov_seq <- allcov_df %>% filter(n_seqs > 0) %>% nrow # n coronaviruses with at least 1 sequence
n_cov_seq_eid2 <- allcov_df %>% filter(n_seqs > 0 & Hostspp > 0) %>% nrow # n coronaviruses with at least 1 sequence and at least 1 EID2 host
n_seq <- allcov_meta_df %>% nrow # n nucleotide sequence entries
n_cds <- cov_enc_df %>% nrow # n coding sequences
# 
# ###########################################################
# # Prep files to run VirusFeatures (R Orton, Java program) #
# ###########################################################
# 
# # Write accesions/taxids/names to file for VirusFeatures
# # SCOPE HERE TO WRITE THE ACCESSION NUMBERS IN EARLIER INTO TAX_MATCHING_DF AND PROCESS THESE WITHOUT EXTRACTING HERE?
# accessions_table <- data.frame(
#   genaccession = unlist(lapply(RefSeq_found_summaries, function(x) extract_from_esummary(x, "oslt")))[c(FALSE,TRUE)],
#   taxid = unlist(lapply(RefSeq_found_summaries, function(x) extract_from_esummary(x, "taxid"))),
#   title = unlist(lapply(RefSeq_found_summaries, function(x) extract_from_esummary(x, "title")))
# )
# 
# # Cross-reference and append taxonomic family of each virus
# accessions_table <- merge(accessions_table, 
#                           tax_matching_df %>% 
#                             select(chosen_refseq_taxid, family_name),
#                           by.x = "taxid",
#                           by.y = "chosen_refseq_taxid",
#                           all.x = TRUE)
# 
# # VirusFeatures program only considers sequences to belong to the same virus if they fall under the same accession number
# # Manually editing accession numbers from here such that sequences are correctly grouped
# 
# # Identify the first accession number for each taxid and add as a new column
# accessions_table <- merge(x=accessions_table, 
#                           y=accessions_table %>% 
#                             select(genaccessionref = genaccession, taxid, -title) %>%
#                             group_by(taxid) %>% 
#                             arrange(genaccessionref) %>% 
#                             slice(1) %>% 
#                             ungroup(),
#                           by="taxid")
# 
# # Read FASTA file back to adjust format
# f <- read.fasta("RefSeq_EID2.fasta")
# 
# # Manually replace existing accession numbers in the fasta
# pattern <- as.character(accessions_table$genaccession)
# replace <- as.character(accessions_table$genaccessionref)
# names(replace) <- pattern
# names(f) <- str_replace_all(names(f), replace)
# # seqinr's write feature can save as single-line FASTA
# write.fasta(f, names(f), nbchar=1e+16, file="VirusFeatures\\RefSeq_EID2_VirusFeatures.fasta")
# 
# # Write a cleaned up version of the accession table
# write.table(subset(accessions_table, genaccession %in% genaccessionref,
#                    select = c(genaccession, taxid, title, family_name)),
#             file="VirusFeatures\\RefSeq_EID2_accessions.txt", 
#             sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)