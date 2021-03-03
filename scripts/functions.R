####################
# Define functions #
####################

# Set up accessible colour blindness palette
cbbPalette_ordered_bw <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#b9b9b9")

# Search GenBank sequences
Seq_searcher <- function(x) {
  Seq_result <- entrez_search(
    db = "nuccore", term = paste0("txid", x, "[Organism:noexp] AND ", searchterm),
    retmax = 10000
  )
}

# Extract Genbank sequence metadata
Seq_summary <- function(x) {

  # Break apart queries to run sequentially so as not to overload API
  query_index <- split(seq(1, length(x)), ceiling(seq_along(seq(1, length(x))) / 300))
  Seq_result <- vector("list", length(query_index))

  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_summary(db = "nuccore", id = x[unlist(query_index[i])])
    Sys.sleep(5) # Delay
  }

  if (length(x) == 1) {
    return(Seq_result)
  } else {
    return(Seq_result %>% flatten() %>% unname())
  }
}

# Extract Genbank sequences as FASTA format
Seq_FASTA <- function(x) {

  # Break apart queries to run sequentially so as not to overload API
  query_index <- split(seq(1, length(x)), ceiling(seq_along(seq(1, length(x))) / 300))
  Seq_result <- vector("list", length(query_index))

  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_fetch(db = "nuccore", id = x[unlist(query_index[i])], rettype = "fasta_cds_na")
    Sys.sleep(5) # Delay
  }
  return(flatten_chr(Seq_result) %>% paste(collapse = ""))
}

# Clean "title" field of GenBank entires to separate out gene/protein names
metadata_title_cleaner <- function(x) {
  x %>%
    as.character() %>%
    strsplit(., "(?<=.)(?= cds)", perl = TRUE) %>% # split title into chunks based on "cds" separator
    lapply(function(x) {
      x[grepl("spike|surface gly|s gly|s prot|S gene|peplom|(S1)| S1 |(S2)| S2 |subunit 1|subunit 2|1 subunit|2 subunit", x, ignore.case = TRUE)] %>% # search for the chunk describing the S protein
        word(., -1)
    }) # save whether partial or complete (which should be last word in string)
}

# Calculate all genome composition bias features
calc_composition_bias <- function(df) {

  # Calculate total nucleotides
  df %<>% mutate(length = A + C + G + T)

  # Calculate amino acid frequencies
  for (i in 1:length(unique(codon_ref$aminoacid))) {
    df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")] <-
      df %>%
      select(subset(codon_ref, aminoacid == unique(codon_ref$aminoacid)[i])$codon) %>%
      rowSums()
  }

  # Calculate total dinucleotides (pos 1-2, pos 2-3, pos 3-1), total codons
  df %<>% mutate(
    n_dinucs = (select(., matches("^[A|C|G|T][A|C|G|T]$")) %>% rowSums()),
    n_dinucs_p1 = (select(., matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% rowSums()),
    n_dinucs_p2 = (select(., matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% rowSums()),
    n_dinucs_p3 = (select(., matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% rowSums()),
    n_codons = (select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% rowSums())
  )

  # Calculate nucleotide biases
  df %<>% mutate_at(vars(matches("^[A|C|G|T]$")), .funs = list(Bias = ~ . / length))

  # Calculate dinucleotide biases
  for (i in 1:ncol(df)) {
    focalcol <- colnames(df)[i]

    if (grepl("^[A|C|G|T][A|C|G|T]$", focalcol)) {
      df[, paste0(focalcol, "_Bias")] <-
        (df[, i] / df$n_dinucs) /
          (df[, substr(focalcol, 1, 1)] / df$length * df[, substr(focalcol, 2, 2)] / df$length)
    }
  }

  # Calculate dinucleotide biases separately for positions 1-2, 2-3, 3-1
  # (Need to calculate nucleotide frequencies for these positions first)

  for (nuc in c("A", "C", "G", "T")) {
    df[paste0(nuc, "_p1")] <- apply(
      df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")), # Select only columns describing dinucleotides at position 1-2
      1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% names() %>% str_count(nuc))
    ) # Use dot product to calculate number of respective nucleotides at position 1-2

    df[paste0(nuc, "_p2")] <- apply(
      df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")), # Select only columns describing dinucleotides at position 2-3
      1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% names() %>% str_count(nuc))
    ) # Use dot product to calculate number of respective nucleotides at position 2-3

    df[paste0(nuc, "_p3")] <- apply(
      df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")), # Select only columns describing dinucleotides at position 3-1
      1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% names() %>% str_count(nuc))
    ) # Use dot product to calculate number of respective nucleotides at position 3-1
  }

  df$n_nucs_p1 <- df %>%
    select(matches("^[A|C|G|T]_p1$")) %>%
    rowSums()
  df$n_nucs_p2 <- df %>%
    select(matches("^[A|C|G|T]_p2$")) %>%
    rowSums()
  df$n_nucs_p3 <- df %>%
    select(matches("^[A|C|G|T]_p3$")) %>%
    rowSums()

  for (i in 1:ncol(df)) {
    focalcol <- colnames(df)[i]

    if (grepl("^[A|C|G|T][A|C|G|T]_p1$", focalcol)) {
      df[, paste0(focalcol, "_Bias")] <-
        (df[, i] / df$n_dinucs_p1) /
          (df[, paste0(substr(focalcol, 1, 1), "_p1")] / df$n_nucs_p1 * df[, paste0(substr(focalcol, 2, 2), "_p1")] / df$n_nucs_p1)
    } else if (grepl("^[A|C|G|T][A|C|G|T]_p2$", focalcol)) {
      df[, paste0(focalcol, "_Bias")] <-
        (df[, i] / df$n_dinucs_p2) /
          (df[, paste0(substr(focalcol, 1, 1), "_p2")] / df$n_nucs_p2 * df[, paste0(substr(focalcol, 2, 2), "_p2")] / df$n_nucs_p2)
    } else if (grepl("^[A|C|G|T][A|C|G|T]_p3$", focalcol)) {
      df[, paste0(focalcol, "_Bias")] <-
        (df[, i] / df$n_dinucs_p3) /
          (df[, paste0(substr(focalcol, 1, 1), "_p3")] / df$n_nucs_p3 * df[, paste0(substr(focalcol, 2, 2), "_p3")] / df$n_nucs_p3)
    }
  }

  # Calculate Relative Synonymous Codon Usage
  for (i in 1:ncol(df)) {
    focalcol <- colnames(df)[i]

    if (grepl("^[A|C|G|T][A|C|G|T][A|C|G|T]$", focalcol)) {
      df[, paste0(focalcol, "_Bias")] <-
        df[, i] * subset(codon_ref, codon == focalcol)$deg /
          (df %>%
            select(subset(codon_ref, aminoacid == subset(codon_ref, codon == focalcol)$aminoacid)$codon) %>%
            rowSums())
    }
  }

  return(df)
}

# Calculate F1micro and F1macro scores
# Sourced from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/ (Matthias Döring, 2018)
get.conf.stats <- function(cm) {
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive]
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positive, x$positive])
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}

get.micro.f1 <- function(cm) {
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp)
  fn <- sum(cm.summary$fn)
  fp <- sum(cm.summary$fp)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  return(f1)
}

get.macro.f1 <- function(cm) {
  c <- cm[[1]]$byClass # a single matrix is sufficient
  re <- sum(c[, "Recall"]) / nrow(c)
  pr <- sum(c[, "Precision"]) / nrow(c)
  f1 <- 2 * ((re * pr) / (re + pr))
  return(f1)
}

# Convenience function for setting up data to plot species-labelled prediction figure
rearrange_to_plot_pred_fig_spp <- function(host) {
  order <- valid_df_raw %>%
    filter(host_category == host) %>%
    group_by(childtaxa_name) %>%
    summarise(correct = mean(!!sym(host))) %>%
    arrange(-correct) %>%
    pull(childtaxa_name)
  valid_df_raw %>%
    filter(host_category == host) %>%
    arrange(factor(childtaxa_name, levels = order), -!!sym(host)) %>%
    return()
}
