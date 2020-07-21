##################
# Setup packages #
##################

rm(list=ls())

# Requires taxizedb_0.1.9.9130 or greater which > CRAN version
#devtools::install_github("ropensci/taxizedb")

library(Biostrings)
library(coRdon)
library(ggbiplot)
library(ggmosaic)
library(gplots)
library(kableExtra)
library(knitr)
library(lattice)
library(magrittr)
library(patchwork)
library(pbapply)
library(plyr)
library(plotly)
library(rentrez)
library(reshape2)
library(rmarkdown)
library(R.utils)
library(seqinr)
library(stringr)
library(taxizedb)
library(tidyverse)
library(VennDiagram)
library(XML)

###############
# Set options #
###############

# Set global variables
setwd("C:\\Users\\Liam\\Desktop\\CoV Genomics")
set_entrez_key("7e8a75a8d92089428e7489e7dc2ea85b0708")

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- data.frame(aminoacid = Biostrings::GENETIC_CODE) %>%
  rownames_to_column("codon") %>%
  mutate(aminoacid = gsub("\\*","X",aminoacid)) %>%     # replace stop codon symbol "*" as "X"
  group_by(aminoacid) %>%
  mutate(deg = n())

# Download and connect to full NCBI taxonomy
db_download_ncbi() # last downloaded to update - 31/3/20
src_ncbi <- src_ncbi()

# Do you want to load previously downloaded data instead of re-extracting?
load_prev_seqs <- TRUE

# Set searchterms
searchterm <- '(spike[Title] OR "S gene"[Title] OR "S protein"[Title] OR "S glycoprotein"[Title] OR "S1 gene"[Title] OR "S1 protein"[Title] OR "S1 glycoprotein"[Title] OR peplomer[Title] OR peplomeric[Title] OR peplomers[All Title] OR "complete genome"[Title]) NOT (patent[Title] OR vaccine OR artificial OR construct OR recombinant[Title])'

# add "surface protein" "surface glycoprotein"?

# Set length filters for whole genome sequences and spike protein sequences
min_wgs_length <- 20000
max_wgs_length <- 32000
min_spike_length <- 2000

# Do you want to truncate sequences that are double-counted through frameshifts? (TRUE = use sequence as-is, FALSE = use sequence as-transcribed)
frameshift_correct <- TRUE

###############
# Run scripts #
###############

header(verbose, "Loading custom functions", padding=0)
source("scripts\\functions.R" )

header(verbose, "Extracting and processing sequence data", padding=0)
source("scripts\\process_cov_seq.R" )

# Save data needed for ML
save(allcov_df, cov_spikes_df, cov_wg_df, file = paste0("cov_ML_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
#save(allcov_df, cov_wg_df, file = paste0("cov_ML_dfs_noframeshift_", format(Sys.time(), "%d_%m_%y"), ".RData"))

# Render lab books
render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.Rmd",
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.html")

# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.html")