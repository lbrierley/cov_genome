##################
# Setup packages #
##################

rm(list=ls())

# Requires taxizedb_0.1.9.9130 or greater which > CRAN version
#devtools::install_github("ropensci/taxizedb")

library(coRdon)
library(ggalluvial)
library(ggmosaic)
library(knitr)
library(lattice)
library(magrittr)
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
setwd("C:\\Users\\Liam\\Desktop\\CoV Genomics\\data")
set_entrez_key("7e8a75a8d92089428e7489e7dc2ea85b0708")
# LOAD PRE SAVED DATA load("TEMP_03_02_20.RData")

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- data.frame(aminoacid = Biostrings::GENETIC_CODE) %>%
  rownames_to_column("codon") %>%
  group_by(aminoacid) %>%
  mutate(deg = n())

# Download and connect to full NCBI taxonomy
db_download_ncbi() # last downloaded to update - 31/3/20
src_ncbi <- src_ncbi()

# Do you want to load previously downloaded data instead of re-extracting?
load_prev_seqs <- TRUE

####################
# Define functions #
####################

# add "surface protein" "surface glycoprotein"?

Seq_searcher <- function(x){
  Seq_result <- entrez_search(db="nuccore", term=paste0('txid', x, '[Organism:noexp] AND (spike[Title] OR "S gene"[Title] OR "S protein"[Title] OR "S glycoprotein"[Title] OR "S1 gene"[Title] OR "S1 protein"[Title] OR "S1 glycoprotein"[Title] OR peplomer[Title] OR peplomeric[Title] OR peplomers[All Title] OR "complete genome"[Title]) NOT (patent[Title] OR vaccine OR artificial OR construct OR recombinant[Title])'),
                              retmax=10000)
}

Seq_summary <- function(x){
  
  query_index <- split(seq(1,length(x)), ceiling(seq_along(seq(1,length(x)))/300))
  Seq_result <- vector("list", length(query_index))
  
  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_summary(db = "nuccore",id = x[unlist(query_index[i])])
    Sys.sleep(5)
  }
  return(flatten(Seq_result))
}

Seq_FASTA <- function(x){
  
  query_index <- split(seq(1,length(x)), ceiling(seq_along(seq(1,length(x)))/300))
  Seq_result <- vector("list", length(query_index))
  
  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_fetch(db = "nuccore",id = x[unlist(query_index[i])], rettype="fasta_cds_na")
    Sys.sleep(5)
  }
  return(flatten_chr(Seq_result) %>% paste(collapse=""))
}