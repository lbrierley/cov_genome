####################
# Define functions #
####################

Seq_searcher <- function(x){
  Seq_result <- entrez_search(db="nuccore", term=paste0('txid', x, '[Organism:noexp] AND ', searchterm),
                              retmax=10000)
}

Seq_summary <- function(x){
  
  query_index <- split(seq(1,length(x)), ceiling(seq_along(seq(1,length(x)))/300))
  Seq_result <- vector("list", length(query_index))
  
  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_summary(db = "nuccore",id = x[unlist(query_index[i])])
    Sys.sleep(5)
  }
  
  if(length(x) == 1){
    return(Seq_result)
  } else {
    return(Seq_result %>% flatten %>% unname)
  }
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

metadata_title_cleaner <- function(x){
  x %>%
    as.character %>%
    strsplit(., "(?<=.)(?= cds)", perl=TRUE) %>%                # split title into chunks based on "cds" separator
    lapply(function(x) x[grepl("spike|surface gly|s gly|s prot|S gene|peplom|(S1)| S1 |(S2)| S2 |subunit 1|subunit 2|1 subunit|2 subunit", x, ignore.case = TRUE)] %>%       # search for the chunk describing the S protein
             word(., -1))                                       # save whether partial or complete (which should be last word in string)
}

genomic_pca <- function(df, vars, outcome, choices = 1:2){
  
  df_name <- deparse(substitute(df))
  
  df %<>% as.data.frame
  
  # Merge in relevant outcome column if it doesn't already exist
  if (!(outcome %in% names(df))){
  df %<>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome)),
                 by = c("taxid" = "childtaxa_id"))
  }
  
  # Create relevant principal components analysis
  
  if (vars == "dinucs"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T]_Bias$")) %>% prcomp
  } else if (vars == "codons"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% prcomp
  } else if (vars == "codons_nostop"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias)) %>% prcomp
  } else if (vars == "aa"){
    pca <- df %>% select(matches("^.*_aa_Bias$")) %>% prcomp
  } else {
    stop("no valid variables selected")
  }
  
  # Write summary
  sink(paste0("figs\\pcasumm_",df_name,"_",vars,".txt"))
  summary(pca)
  sink()
  
  # Screeplot
  ggscreeplot(pca) +
    geom_hline(yintercept=1/length(pca$sdev), alpha = 0.5, color="dodgerblue", lty="dashed", size=1.5) +
    theme_bw() +
    ggsave(paste0("figs\\scree_",df_name,"_",vars,".png"), width = 8, height = 5)
  
  # Biplot
  g <- ggbiplot(pca,
                choices = choices,
                groups = df[, outcome],
                ellipse = TRUE,
                alpha = 0.4,
                varname.abbrev=TRUE) +
    geom_point(alpha=0, aes(fill= df[, outcome], label=df$childtaxa_name)) +
    theme(legend.position='none') +
    theme_bw()
  ggplotly(g) %>% hide_legend()
  
}


calc_composition_bias <- function(df){
  
  # Calculate total nucleotides
  df %<>% mutate(length = A+C+G+T)
  
  # Calculate amino acid frequencies
  for (i in 1:length(unique(codon_ref$aminoacid))) {
    df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")] <- 
      df %>% select(subset(codon_ref, aminoacid == unique(codon_ref$aminoacid)[i])$codon) %>% rowSums()
  }
  
  # Calculate total dinucleotides (pos 1-2, pos 2-3, pos 3-1), total codons
  df %<>% mutate(n_dinucs = (select(., matches("^[A|C|G|T][A|C|G|T]$")) %>% rowSums),
                 n_dinucs_p1 = (select(., matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% rowSums),
                 n_dinucs_p2 = (select(., matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% rowSums),
                 n_dinucs_p3 = (select(., matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% rowSums),
                 n_codons = (select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% rowSums))
  
  # Calculate nucleotide biases
  df %<>% mutate_at(vars(matches("^[A|C|G|T]$")), .funs = list(Bias = ~./length))
  
  
  # # MUST BE A DPLYR MUTATE AT WAY OF DOING THIS
  # cocoputs_dinucs %>%
  #   mutate_at(vars(matches("^[A|T|G|C|U]p[A|T|G|C|U]$")), funs(
  #     (./X..Dinucleotides)/ # Numerator, Nxy/Dtot
  #       (!!parse_quosure(deparse(substitute(.)) %>% substr(1,1))/X..Nucleotides * # Denominator, Nx/Ntot, extracting nucleotide X from col name
  #          !!parse_quosure(deparse(substitute(.)) %>% substr(3,3))/X..Nucleotides) # Denominator, Ny/Ntot, extracting nucleotide Y from col name
  #   ))
  # # Should work but object v not found? Comes back to nonstandard evaluation, a topic that seems tricky to get a full handle on
  
  # Calculate dinucleotide biases
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T]$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs)/
        (df[,substr(focalcol,1,1)]/df$length * df[,substr(focalcol,2,2)] / df$length)
    }
  }
  
  # Calculate dinucleotide biases separately for positions 1-2, 2-3, 3-1
  # But first, need to calculate nucleotide frequencies for these positions - slow but easiest way of doing the calculation
  
  for (nuc in c("A","C","G","T")){
    df[paste0(nuc,"_p1")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")),     # Select only columns describing dinucleotides at position 1-2
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 1-2
    
    df[paste0(nuc,"_p2")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")),     # Select only columns describing dinucleotides at position 2-3
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 2-3
    
    df[paste0(nuc,"_p3")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")),     # Select only columns describing dinucleotides at position 3-1
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 3-1
  }
  
  df$n_nucs_p1 <- df %>% select(matches("^[A|C|G|T]_p1$")) %>% rowSums
  df$n_nucs_p2 <- df %>% select(matches("^[A|C|G|T]_p2$")) %>% rowSums
  df$n_nucs_p3 <- df %>% select(matches("^[A|C|G|T]_p3$")) %>% rowSums
  
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T]_p1$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p1)/
        (df[,paste0(substr(focalcol,1,1),"_p1")]/df$n_nucs_p1 * df[,paste0(substr(focalcol,2,2),"_p1")] / df$n_nucs_p1)
    } else if (grepl("^[A|C|G|T][A|C|G|T]_p2$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p2)/
        (df[,paste0(substr(focalcol,1,1),"_p2")]/df$n_nucs_p2 * df[,paste0(substr(focalcol,2,2),"_p2")] / df$n_nucs_p2)
    } else   if (grepl("^[A|C|G|T][A|C|G|T]_p3$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p3)/
        (df[,paste0(substr(focalcol,1,1),"_p3")]/df$n_nucs_p3 * df[,paste0(substr(focalcol,2,2),"_p3")] / df$n_nucs_p3)
    }
  }
  
  # Calculate Relative Synonymous Codon Usage
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T][A|C|G|T]$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        df[,i]*subset(codon_ref, codon == focalcol)$deg/
        (df %>%
           select(subset(codon_ref, aminoacid == subset(codon_ref, codon == focalcol)$aminoacid)$codon) %>%
           rowSums())
    }
  }
  
  
  # Calculate amino acid biases - denominator uses total amino acids, including stop codons, so can just use total codons
  for (i in 1:length(unique(codon_ref$aminoacid))) {
    df[, paste0(unique(codon_ref$aminoacid)[i], "_aa_Bias")] <- 
      df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")]/df$n_codons
  }
  
  
  # # Calc codon pair bias as log (freq codon pair)/((freq codon A*freq codon B/freq aacid A*freq aacid B) * freq aacid pair) following Coleman et al. 2008
  # # WORKS - BUT NOT CURRENTLY USING
  #
  # # Calculate codon pair biases
  # for (i in 1:nrow(codon_ref)) {
  #   for (j in 1:nrow(codon_ref)) {
  #     
  #     # Calculate frequency of pairs of corresponding amino acids for codons i and j  
  #     aminoacidpaircounts <- df %>% 
  #       select(do.call(paste0, 
  #                      expand.grid(codon_ref %>% subset(aminoacid == codon_ref$aminoacid[i]) %>% .$codon,
  #                                  codon_ref %>% subset(aminoacid == codon_ref$aminoacid[j]) %>% .$codon))) %>% rowSums()
  #     
  #     df[, paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <- 
  #       log(
  #         df[, paste0(codon_ref$codon[i], codon_ref$codon[j])]/
  #           (((df[, codon_ref$codon[i]]*df[, codon_ref$codon[j]])/
  #               (df[, paste0(codon_ref$aminoacid[i], "_aa")]*df[, paste0(codon_ref$aminoacid[j], "_aa")]))*
  #              aminoacidpaircounts))
  #     
  #     # If frequency of codon pair = 0 but frequency of amino acid pair != 0 specifiy codon pair bias as NA as a marker to replace later
  #     df[which(aminoacidpaircounts != 0 & df[, paste0(codon_ref$codon[i], codon_ref$codon[j])] == 0), paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <- NA
  #     
  #   }
  # }
  # 
  # # Below calculations are easily changeable!
  # # Work out column of mean codon pair bias per virus across all non-NaN and non-NA values
  # df %<>% mutate(mean_CPB = rowMeans(select(., (matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"))), na.rm=TRUE))
  # 
  # # Do it not including pairs involving stop codons
  # df %<>% mutate(mean_CPB_nostop = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA|TAG|TAA")), na.rm=TRUE))
  # df %<>% mutate(mean_CPB_nostop1 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("^TGA|^TAG|^TAA")), na.rm=TRUE))
  # df %<>% mutate(mean_CPB_nostop2 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA.Bias|TAG.Bias|TAA.Bias")), na.rm=TRUE))
  # 
  # # Following Babayan et al. rules: if frequency of amino acid pair = 0 replace NaN with mean codon pair bias
  # df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")), ~ifelse(is.nan(.), mean_CPB, .))
  # 
  # # Following Babayan et al. rules: if frequency of codon pair = 0 but frequency of amino acid pair != 0 replace NA with -9999 to indicate extreme underrepresentation
  # df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")),~ifelse(is.na(.), -9999, .))
  
  return(df)
}