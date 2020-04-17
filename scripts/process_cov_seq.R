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

# Append genus name, sequence counts of viruses found
allcov_df <- data.frame(allcov_df, 
                        n_seqs = lapply(allcov_results, function(x) x$count) %>% unlist())

genus_list <- lapply(allcov_df$childtaxa_id %>% classification, function(x) x %>% filter(rank == "genus") %>% .$name)
is.na(genus_list) <- lengths(genus_list) == 0

allcov_df$genus <- genus_list %>% unlist %>% unname

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

# Specify SARS-CoV-2 as a human virus, which currently isn't recorded in EID2
allcov_df[allcov_df$childtaxa_id == 2697049,"h_human"] <- "1"
allcov_df[allcov_df$childtaxa_id == 2697049,"n_groups"] <- 1
allcov_df[allcov_df$childtaxa_id == 2697049,"Hostspp"] <- 1

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

# Read FASTA back in in coRdon format and calculate codon table
cov_cord <- readSet(file="data\\allcov_GenBank.fasta")

# Construct summary df for individual cds
cov_enc_df <- data.frame(title = cov_cord %>% names,
                         accessionversion = cov_cord %>% names %>% str_match("lcl\\|(.*?)_cds_") %>% .[,2] %>% as.character,
                         gene = cov_cord %>% names %>% str_match("\\[gene=(.*?)\\]") %>% .[,2],
                         protein = cov_cord %>% names %>% str_match("\\[protein=(.*?)\\]") %>% .[,2],
                         length = cov_cord %>% width,
                         enc = cov_cord %>% codonTable %>% ENC(stop.rm=FALSE),              # Calculate Effective Number of Codons (including STOP codons)
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
) %>% rename_at(vars(G.C), ~ "GC_content")

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

cov_enc_df %<>% mutate(include = case_when(grepl("trunc|mutant", protein, ignore.case=TRUE) ~ "N")) %>%
  replace_na(list(include = "Y"))

## Code to assist choosing regexp for assigning seqtype
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 0) %>% kable
# cov_enc_df %>% filter(is.na(seqtype)) %>% select(protein) %>% table %>% as.data.frame %>% filter(Freq > 100) %>% arrange(Freq) %>% kable #most common protein labels not captured
# allcov_meta_df %>% filter(as.character(accessionversion) %in% (cov_enc_df %>% filter(is.na(protein)) %>% .$accession %>% unique)) %>% nrow # number of covs with unlabelled spike proteins

# Add in taxids corresponding to accesion nos - plyr::join is used to retain row order
cov_enc_df <- plyr::join(allcov_meta_df %>% select(accessionversion, taxid, complete),
                    cov_enc_df,
                    by = "accessionversion",
                    all.y = TRUE)

# Assign partial if labelled as complete spike but too small to be complete spike
cov_enc_df %<>% mutate(complete = ifelse(seqtype == "S" & length < 2000, "partial_spike", complete))

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

# Calculate mean length of whole genomes per coronavirus taxid and merge with main dataset

allcov_df <- merge(allcov_df,
                   allcov_meta_df %>% 
                     filter(complete == "whole_genome") %>% 
                     group_by(taxid) %>% 
                     summarise(mean_wgs_length = mean(length, na.rm=TRUE)) %>% 
                     as.data.frame,
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

# Calculate mean enc, GC content of complete genomes per coronavirus taxid and merge with main dataset

wg_df <- cov_enc_df %>% 
  filter(include == "Y") %>% 
  group_by(accessionversion) %>% 
  summarise_at(vars(matches("^[A|C|G|T]$|^[A|C|G|T][A|C|G|T]$|^[A|C|G|T][A|C|G|T][A|C|G|T]$")), funs(sum)) %>% # Calculate nuc, dinuc, codon counts over entire sequence (including all proteins)
  mutate(GC_content = 100*(G+C)/(A+C+G+T))

wg_df$enc <- wg_df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% codonTable() %>% ENC

allcov_df <- merge(allcov_df,
                   merge(cov_enc_df %>% select(taxid, accessionversion),
                         wg_df,
                         by = "accessionversion",
                         all.x = TRUE) %>% 
                     filter(accessionversion %in% (allcov_meta_df %>% filter(complete == "whole_genome") %>% .$accessionversion)) %>%
                     group_by(taxid) %>% 
                     summarise(mean_wgs_enc = mean(enc, na.rm=TRUE),
                               mean_wgs_GC = mean(GC_content, na.rm=TRUE)),
                   by.x = "childtaxa_id",
                   by.y = "taxid",
                   all.x = TRUE)

# Better to plot length/GC/ENC of ALL whole genomes than mean within taxid?



# ### why is mean_wgs_enc != mean_wgs_length in terms of NAs? it's those missing sequences again - probably want to investigate this a little as it is 78 viruses!
# hi <- allcov_meta_df %>% filter(complete == "whole_genome") %>% distinct(accessionversion) %>% .$accessionversion
# # 2954
# lo <- wg_df %>% filter(accessionversion %in% (allcov_meta_df %>% filter(complete == "whole_genome") %>% .$accessionversion)) %>% .$accessionversion
# # 2876
# allcov_meta_df %>% filter(accessionversion %in% hi[!(hi %in% lo)]) %>% .$uid %>%  entrez_fetch(db = "nuccore",id = ., rettype="fasta_cds_na")
# 
# # length is pulled from metadata,  whereas GC/ENC must be taken from sequence, so even if it returns "/n" we can still get length
# allcov_df %>% filter(!is.na(mean_wgs_length)) %>% nrow
# allcov_df %>% filter(!is.na(mean_wgs_GC)) %>% nrow
# 
# missingids <- allcov_df %>% filter(!is.na(mean_wgs_length) & is.na(mean_wgs_GC)) %>% .$childtaxa_id
# missingaccession <- allcov_meta_df %>% filter(taxid %in% missingids & complete == "whole_genome") %>% .$accessionversion
# allcov_meta_df %>% filter(accessionversion %in% missingaccession) %>% .$uid
# 
# # just pulling in "\n"...  THE REASON IS THE SEQUENCES ARE NOT ACTUALLY ANNOTATED WITH CODING SEQUENCE REGIONS AND NONCODING SEQUENCE REGIONS
# IF WE REALLY NEED THESE I SUPPOSE WE COULD TRY ALIGNING?
# CAN THEY BE FILTERED EARLIER TO KEEP NUMBERS UPDATED, E.G. IN ENTREZ_SEARCH AND ENTREZ_SUMMARY

#################################
# Calculate values of interest  #
#################################

# Count total sequences for various criteria
n_cov <- nrow(allcov_df) # n coronaviruses
n_cov_seq <- allcov_df %>% filter(n_seqs > 0) %>% nrow # n coronaviruses with at least 1 sequence
n_cov_seq_eid2 <- allcov_df %>% filter(n_seqs > 0 & Hostspp > 0) %>% nrow # n coronaviruses with at least 1 sequence and at least 1 EID2 host
n_seq <- allcov_meta_df %>% nrow # n nucleotide sequence entries
n_cds <- cov_enc_df %>% nrow # n coding sequences

# # 492 viruses should have complete spike protein based on metadata, but only 420 do based on the cds labels in the FASTA!miss_spike <- tab492[!(tab492 %in% tab420)]
# tab420 <- allcov_df %>% filter(!is.na(n_S) & n_S > 0) %>% .$childtaxa_id
# tab492 <-  allcov_meta_df %>% group_by(taxid, complete) %>% count(complete) %>% spread(complete, n, fill=0) %>% rename_at(vars(-taxid), ~ paste0("n_",.)) %>% filter(n_complete_spike > 0 | n_whole_genome > 0) %>% .$taxid %>% as.character
# 
# miss_spike <- tab492[!(tab492 %in% tab420)] # Identify taxids missing the spike protein
# 
# allcov_df %>% filter(childtaxa_id %in% miss_spike) %>% .$n_other %>% table(exclude = FALSE) # 8 come through in other sequences, but 64 have NAs
# 
# # For the 8, metadata describes them as complete genome...
# 
# allcov_meta_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & !(is.na(n_other))) %>% .$childtaxa_id)) %>% .$title
# 
# # ...yet they don't have any spike proteins (or none annotated)!
# 
# cov_enc_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & !(is.na(n_other))) %>% .$childtaxa_id)) %>% .$protein
# 
# # For the 64 - entrez_fetch just.. fails to fetch them and returns blanks???
# 
# allcov_meta_df %>% filter(taxid %in% (allcov_df %>% filter(childtaxa_id %in% miss_spike & is.na(n_other)) %>% .$childtaxa_id)) %>% .$uid %>% as.character %>% entrez_fetch(db = "nuccore",id = ., rettype="fasta_cds_na")
# 
# # So 420 is the correct number of viruses with a spike protein...
# sum(allcov_df$n_S, na.rm=TRUE) # 5213 individual spike protein sequences

###################################
# Establish spike protein dataset #
###################################

# Filter dataset to final complete spike protein dataset
cov_spikes_df <- cov_enc_df %>% filter(seqtype == "S" & include == "Y" & complete != "partial_spike")

cov_enc_df[which(cov_enc_df$seqtype == "S" & cov_enc_df$include == "Y" & cov_enc_df$complete != "partial_spike"), ]

# Write FASTA of final complete spike protein dataset
cov_cord 40801

#######################################
# Calculate genome composition biases #
#######################################

# Calculate amino acid frequencies
for (i in 1:length(unique(codon_ref$aminoacid))) {
  cov_spikes_df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")] <- 
    cov_spikes_df %>% select(subset(codon_ref, aminoacid == unique(codon_ref$aminoacid)[i])$codon) %>% rowSums()
}

# Calculate total dinucleotides (pos 1-2, pos 2-3, pos 3-1), total codons
cov_spikes_df %<>% mutate(n_dinucs = (select(., matches("^[A|C|G|T][A|C|G|T]$")) %>% rowSums),
                          n_dinucs_p1 = (select(., matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% rowSums),
                          n_dinucs_p2 = (select(., matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% rowSums),
                          n_dinucs_p3 = (select(., matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% rowSums),
                          n_codons = (select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% rowSums))

# Calculate nucleotide biases
cov_spikes_df %<>% mutate_at(vars(matches("^[A|C|G|T]$")), .funs = list(Bias = ~./length))


# # MUST BE A DPLYR MUTATE AT WAY OF DOING THIS
# cocoputs_dinucs %>%
#   mutate_at(vars(matches("^[A|T|G|C|U]p[A|T|G|C|U]$")), funs(
#     (./X..Dinucleotides)/ # Numerator, Nxy/Dtot
#       (!!parse_quosure(deparse(substitute(.)) %>% substr(1,1))/X..Nucleotides * # Denominator, Nx/Ntot, extracting nucleotide X from col name
#          !!parse_quosure(deparse(substitute(.)) %>% substr(3,3))/X..Nucleotides) # Denominator, Ny/Ntot, extracting nucleotide Y from col name
#   ))
# # Should work but object v not found? Comes back to nonstandard evaluation, a topic that seems tricky to get a full handle on

# Calculate dinucleotide biases
for (i in 1:ncol(cov_spikes_df)) {
  
  focalcol <- colnames(cov_spikes_df)[i]
  
  if (grepl("^[A|C|G|T][A|C|G|T]$",focalcol)){
    cov_spikes_df[, paste0(focalcol, "_Bias")] <- 
      (cov_spikes_df[,i]/cov_spikes_df$n_dinucs)/
      (cov_spikes_df[,substr(focalcol,1,1)]/cov_spikes_df$length * cov_spikes_df[,substr(focalcol,2,2)] / cov_spikes_df$length)
  }
}

# Calculate dinucleotide biases separately for positions 1-2, 2-3, 3-1
# But first, need to calculate nucleotide frequencies for these positions
for (nuc in c("A","C","G","T")){
  cov_spikes_df[paste0(nuc,"_p1")] <- apply(cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")),     # Select only columns describing dinucleotides at position 1-2
                                            1, function(x) x %*% (cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 1-2
  
  cov_spikes_df[paste0(nuc,"_p2")] <- apply(cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")),     # Select only columns describing dinucleotides at position 2-3
                                            1, function(x) x %*% (cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 2-3
  
  cov_spikes_df[paste0(nuc,"_p3")] <- apply(cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")),     # Select only columns describing dinucleotides at position 3-1
                                            1, function(x) x %*% (cov_spikes_df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 3-1
}

cov_spikes_df$n_nucs_p1 <- cov_spikes_df %>% select(matches("^[A|C|G|T]_p1$")) %>% rowSums
cov_spikes_df$n_nucs_p2 <- cov_spikes_df %>% select(matches("^[A|C|G|T]_p2$")) %>% rowSums
cov_spikes_df$n_nucs_p3 <- cov_spikes_df %>% select(matches("^[A|C|G|T]_p3$")) %>% rowSums

for (i in 1:ncol(cov_spikes_df)) {
  
  focalcol <- colnames(cov_spikes_df)[i]
  
  if (grepl("^[A|C|G|T][A|C|G|T]_p1$",focalcol)){
    cov_spikes_df[, paste0(focalcol, "_Bias")] <- 
      (cov_spikes_df[,i]/cov_spikes_df$n_dinucs_p1)/
      (cov_spikes_df[,paste0(substr(focalcol,1,1),"_p1")]/cov_spikes_df$n_nucs_p1 * cov_spikes_df[,paste0(substr(focalcol,2,2),"_p1")] / cov_spikes_df$n_nucs_p1)
  } else if (grepl("^[A|C|G|T][A|C|G|T]_p2$",focalcol)){
    cov_spikes_df[, paste0(focalcol, "_Bias")] <- 
      (cov_spikes_df[,i]/cov_spikes_df$n_dinucs_p2)/
      (cov_spikes_df[,paste0(substr(focalcol,1,1),"_p2")]/cov_spikes_df$n_nucs_p2 * cov_spikes_df[,paste0(substr(focalcol,2,2),"_p2")] / cov_spikes_df$n_nucs_p2)
  } else   if (grepl("^[A|C|G|T][A|C|G|T]_p3$",focalcol)){
    cov_spikes_df[, paste0(focalcol, "_Bias")] <- 
      (cov_spikes_df[,i]/cov_spikes_df$n_dinucs_p3)/
      (cov_spikes_df[,paste0(substr(focalcol,1,1),"_p3")]/cov_spikes_df$n_nucs_p3 * cov_spikes_df[,paste0(substr(focalcol,2,2),"_p3")] / cov_spikes_df$n_nucs_p3)
  }
}


# Calculate Relative Synonymous Codon Usage
for (i in 1:ncol(cov_spikes_df)) {
  
  focalcol <- colnames(cov_spikes_df)[i]
  
  if (grepl("^[A|C|G|T][A|C|G|T][A|C|G|T]$",focalcol)){
    cov_spikes_df[, paste0(focalcol, "_Bias")] <- 
      cov_spikes_df[,i]*subset(codon_ref, codon == focalcol)$deg/
      (cov_spikes_df %>%
         select(subset(codon_ref, aminoacid == subset(codon_ref, codon == focalcol)$aminoacid)$codon) %>%
         rowSums())
  }
}


# Calculate amino acid biases - denominator uses total amino acids, including stop codons, so can just use total codons
for (i in 1:length(unique(codon_ref$aminoacid))) {
  cov_spikes_df[, paste0(unique(codon_ref$aminoacid)[i], "_aa_Bias")] <- 
    cov_spikes_df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")]/cov_spikes_df$n_codons
}


# # Calc codon pair bias as log (freq codon pair)/((freq codon A*freq codon B/freq aacid A*freq aacid B) * freq aacid pair) following Coleman et al. 2008
# # WORKS - BUT NOT CURRENTLY USING
#
# # Calculate codon pair biases
# for (i in 1:nrow(codon_ref)) {
#   for (j in 1:nrow(codon_ref)) {
#     
#     # Calculate frequency of pairs of corresponding amino acids for codons i and j  
#     aminoacidpaircounts <- cov_spikes_df %>% 
#       select(do.call(paste0, 
#                      expand.grid(codon_ref %>% subset(aminoacid == codon_ref$aminoacid[i]) %>% .$codon,
#                                  codon_ref %>% subset(aminoacid == codon_ref$aminoacid[j]) %>% .$codon))) %>% rowSums()
#     
#     cov_spikes_df[, paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <- 
#       log(
#         cov_spikes_df[, paste0(codon_ref$codon[i], codon_ref$codon[j])]/
#           (((cov_spikes_df[, codon_ref$codon[i]]*cov_spikes_df[, codon_ref$codon[j]])/
#               (cov_spikes_df[, paste0(codon_ref$aminoacid[i], "_aa")]*cov_spikes_df[, paste0(codon_ref$aminoacid[j], "_aa")]))*
#              aminoacidpaircounts))
#     
#     # If frequency of codon pair = 0 but frequency of amino acid pair != 0 specifiy codon pair bias as NA as a marker to replace later
#     cov_spikes_df[which(aminoacidpaircounts != 0 & cov_spikes_df[, paste0(codon_ref$codon[i], codon_ref$codon[j])] == 0), paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <- NA
#     
#   }
# }
# 
# # Below calculations are easily changeable!
# # Work out column of mean codon pair bias per virus across all non-NaN and non-NA values
# cov_spikes_df %<>% mutate(mean_CPB = rowMeans(select(., (matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"))), na.rm=TRUE))
# 
# # Do it not including pairs involving stop codons
# cov_spikes_df %<>% mutate(mean_CPB_nostop = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA|TAG|TAA")), na.rm=TRUE))
# cov_spikes_df %<>% mutate(mean_CPB_nostop1 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("^TGA|^TAG|^TAA")), na.rm=TRUE))
# cov_spikes_df %<>% mutate(mean_CPB_nostop2 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA.Bias|TAG.Bias|TAA.Bias")), na.rm=TRUE))
# 
# # Following Babayan et al. rules: if frequency of amino acid pair = 0 replace NaN with mean codon pair bias
# cov_spikes_df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")), ~ifelse(is.nan(.), mean_CPB, .))
# 
# # Following Babayan et al. rules: if frequency of codon pair = 0 but frequency of amino acid pair != 0 replace NA with -9999 to indicate extreme underrepresentation
# cov_spikes_df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")),~ifelse(is.na(.), -9999, .))

###########################
# Correspondence analysis #
###########################

# Append genus and virus name information
cov_pca_df <- merge(cov_spikes_df,
                    allcov_df %>% select(childtaxa_id, childtaxa_name, genus, h_human) %>% mutate(genus = gsub("0","unclassified",genus)),
                    all.x = TRUE,
                    by.x = "taxid",
                    by.y = "childtaxa_id")

# Dinucs

dinucs_pca <- cov_pca_df %>% select(matches("^[A|C|G|T][A|C|G|T]_Bias$")) %>% prcomp

summary(dinucs_pca)

ggscreeplot(dinucs_pca) +
  geom_hline(yintercept=1/length(dinucs_pca$sdev), alpha = 0.5, color="dodgerblue", lty="dashed", size=1.5) +
  theme_bw()

# Codons

codons_pca <- cov_pca_df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% prcomp

summary(codons_pca)

ggscreeplot(codons_pca) +
  geom_hline(yintercept=1/length(codons_pca$sdev), alpha = 0.5, color="dodgerblue", lty="dashed", size=1.5) +
  theme_bw()

# Biplots
genome_biplot(codons_pca, "genus")
genome_biplot(codons_pca, "h_human")

# Amino acids

aa_pca <- cov_pca_df %>% select(matches("^.*_aa_Bias$")) %>% prcomp

summary(aa_pca)

ggscreeplot(aa_pca) +
  scale_x_continuous(minor_breaks = seq(0, length(aa_pca$sdev), 5)) +
  geom_hline(yintercept=1/length(aa_pca$sdev), alpha = 0.5, color="dodgerblue", lty="dashed", size=1.5) +
  theme_bw()

# Biplots
genome_biplot(aa_pca, "genus")
genome_biplot(aa_pca, "h_human")
