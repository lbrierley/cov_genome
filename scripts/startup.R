#############################################################
# Supporting R scripts for: Brierley & Fowler et al. 2020   #
#      Predicting the animal hosts of coronaviruses from    #
#   compositional biases of spike protein and whole genome  #
#             sequences through machine learning            #
#                                                           #
#                  Compiled by L. Brierley                  #
#               University of Liverpool, 2021               #
#############################################################


##################
# Setup packages #
##################

rm(list = ls())

# Requires taxizedb_0.1.9.9130 or greater which > CRAN version
# devtools::install_github("ropensci/taxizedb")

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
set_entrez_key("APIKEYGOESHERE") # users should provide their own entrez API key!

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- data.frame(aminoacid = Biostrings::GENETIC_CODE) %>%
  rownames_to_column("codon") %>%
  mutate(aminoacid = gsub("\\*", "X", aminoacid)) %>% # replace stop codon symbol "*" as "X"
  group_by(aminoacid) %>%
  mutate(deg = n())

# Download and connect to full NCBI taxonomy (last downloaded 31/3/20)
db_download_ncbi()
src_ncbi <- src_ncbi()

# Do you want to load previously downloaded data instead of re-extracting?
load_prev_seqs <- TRUE

# Set searchterms
searchterm <-
  '(spike[Title] OR "S gene"[Title] OR "S protein"[Title] OR "S glycoprotein"[Title] OR "S1 gene"[Title] OR "S1 protein"[Title] OR "S1 glycoprotein"[Title] OR peplomer[Title] OR peplomeric[Title] OR peplomers[All Title] OR "complete genome"[Title]) NOT (patent[Title] OR vaccine OR artificial OR construct OR recombinant[Title])'

# Set length filters for whole genome sequences and spike protein sequences
min_wgs_length <- 20000
max_wgs_length <- 32000
min_spike_length <- 2000
min_s1_length <- 1250

# Do you want to truncate sequences that are double-counted through frameshifts?
# (TRUE = use sequence as-is, FALSE = use sequence as-transcribed)
frameshift_correct <- TRUE

# Set options for machine learning analyses
outcome_name <- "group_name"
use_stop_codons <- TRUE
exclude_zoonotic <- TRUE

###############
# Run scripts #
###############

header(verbose, "Loading custom functions", padding = 0)
source("scripts\\functions.R")

header(verbose, "Extracting and processing sequence data", padding = 0)
source("scripts\\process_cov_seq.R")

# header(verbose, "Training random forest models", padding = 0)
# source("scripts\\build_random_forests.R")

# header(verbose, "Calculating partial dependence profiles for most important features", padding = 0)
# source("scripts\\partial_dependence.R")

# header(verbose, "Training random forest models using alternative resampling methods", padding = 0)
# source("scripts\\build_random_forests_samplers.R")

header(verbose, "Constructing figures and tables", padding = 0)
source("scripts\\figures_tables.R")
