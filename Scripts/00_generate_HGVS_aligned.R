library(Biostrings)

# Amino acid 3-letter codes
aa3 <- c(
  A="Ala", R="Arg", N="Asn", D="Asp", C="Cys",
  Q="Gln", E="Glu", G="Gly", H="His", I="Ile",
  L="Leu", K="Lys", M="Met", F="Phe", P="Pro",
  S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val"
)

generate_HGVS_aligned <- function(ref_seq, var_seq,
                                  mutation_type = NULL,
                                  subst_matrix = "BLOSUM62",
                                  gapOpening = -10, gapExtension = -0.5) {
  ref_seq <- toupper(ref_seq)
  var_seq <- toupper(var_seq)
  n_ref <- nchar(ref_seq)
  n_var <- nchar(var_seq)
  
  if (ref_seq == var_seq) return("WT")
  
  ref <- AAString(ref_seq)
  var <- AAString(var_seq)
  
  aln <- pairwiseAlignment(ref, var,
                           substitutionMatrix = subst_matrix,
                           gapOpening = gapOpening, gapExtension = gapExtension,
                           type = "global")
  
  ref_aln <- strsplit(as.character(pattern(aln)), "")[[1]]
  var_aln <- strsplit(as.character(subject(aln)), "")[[1]]
  L <- length(ref_aln)
  
  ## Map reference positions across alignment positions (1..L)
  ref_pos_map <- integer(L)  # ref_pos_map[i] = reference position at alignment column i (0 if gap)
  var_pos_map <- integer(L)
  rp <- 0; vp <- 0
  for (i in seq_len(L)) {
    if (ref_aln[i] != "-") rp <- rp + 1
    if (var_aln[i] != "-") vp <- vp + 1
    ref_pos_map[i] <- rp
    var_pos_map[i] <- vp
  }
  
  # first/last aligned reference positions (for truncation detection)
  idx_var_non_gap <- which(var_aln != "-")
  first_ref_aligned_pos <- if (length(idx_var_non_gap) > 0) ref_pos_map[idx_var_non_gap[1]] else NA_integer_
  last_ref_aligned_pos  <- if (length(idx_var_non_gap) > 0) ref_pos_map[idx_var_non_gap[length(idx_var_non_gap)]] else NA_integer_
  
  # --- Dedicated truncation chunk (activated when mutation_type indicates truncation) ---
  if (!is.null(mutation_type) && grepl("trunc", mutation_type, ignore.case = TRUE)) {
    changes <- character(0)
    is_nterm_req <- grepl("n-?term|nterm|n-terminal|^n\\b", mutation_type, ignore.case = TRUE)
    is_cterm_req <- grepl("c-?term|cterm|c-terminal|^c\\b", mutation_type, ignore.case = TRUE)
    
    # Special case: variant aligns to nothing (completely deleted)
    if (all(var_aln == "-")) {
      changes <- c(changes, paste0("p.", aa3[substr(ref_seq,1,1)], "1_",
                                   aa3[substr(ref_seq,n_ref,n_ref)], n_ref, "del"))
      return(paste(unique(changes), collapse = ";"))
    }
    
    # N-terminal truncation (missing leading residues)
    if (is.na(first_ref_aligned_pos)) {
      # no aligned ref positions -> full deletion handled above
      first_ref_aligned_pos <- NA
    }
    if (is.null(is_nterm_req) || (!is_nterm_req && !is_cterm_req)) {
      # if the user supplied only "truncation" (no N/C), treat as generic: compute both N and C if present
      is_nterm_req <- TRUE
      is_cterm_req <- TRUE
    }
    
    if (is_nterm_req && !is.na(first_ref_aligned_pos) && first_ref_aligned_pos > 1) {
      del_start <- 1
      del_end <- first_ref_aligned_pos - 1
      changes <- c(changes, paste0("p.", aa3[substr(ref_seq,del_start,del_start)], del_start,
                                   "_", aa3[substr(ref_seq,del_end,del_end)], del_end, "del"))
    }
    
    # C-terminal truncation (missing trailing residues)
    if (is_cterm_req && !is.na(last_ref_aligned_pos) && last_ref_aligned_pos < n_ref) {
      del_start <- last_ref_aligned_pos + 1
      del_end <- n_ref
      changes <- c(changes, paste0("p.", aa3[substr(ref_seq,del_start,del_start)], del_start,
                                   "_", aa3[substr(ref_seq,del_end,del_end)], del_end, "del"))
    }
    
    return(paste(unique(changes), collapse = ";"))
  }
  # --- end truncation chunk ---
  
  
  # --- Main path (non-truncation or mutation_type not set to truncation) ---
  changes <- character(0)
  i <- 1
  ref_pos <- 0
  var_pos <- 0
  del_group <- NULL
  del_start <- NULL
  
  while (i <= L) {
    r <- ref_aln[i]
    v <- var_aln[i]
    
    if (r != "-") ref_pos <- ref_pos + 1
    if (v != "-") var_pos <- var_pos + 1
    
    # Substitution
    if (r != "-" && v != "-" && r != v) {
      changes <- c(changes, paste0("p.", aa3[r], ref_pos, aa3[v]))
    }
    
    # Insertion (group contiguous)
    if (r == "-" && v != "-") {
      ins_seq <- c()
      ins_start_pos <- ref_pos # position to the left of insertion (0 means before first residue)
      while (i <= L && ref_aln[i] == "-" && var_aln[i] != "-") {
        ins_seq <- c(ins_seq, var_aln[i])
        i <- i + 1
      }
      ins_3 <- paste(sapply(ins_seq, function(x) aa3[x]), collapse = "")
      if (ins_start_pos == 0) {
        # N-terminal insertion
        first_res <- substr(ref_seq, 1, 1)
        changes <- c(changes, paste0("p.0_", aa3[first_res], "1ins", ins_3))
      } else if (ins_start_pos >= n_ref) {
        # C-terminal insertion (after last residue)
        last_res <- substr(ref_seq, n_ref, n_ref)
        changes <- c(changes, paste0("p.", aa3[last_res], n_ref, "_Ins", ins_3))
      } else {
        right_ref <- min(ins_start_pos + 1, n_ref)
        changes <- c(changes, paste0("p.", aa3[substr(ref_seq, ins_start_pos, ins_start_pos)],
                                     ins_start_pos, "_", aa3[substr(ref_seq, right_ref, right_ref)],
                                     right_ref, "ins", ins_3))
      }
      next
    }
    
    # Deletion (group contiguous)
    if (r != "-" && v == "-") {
      if (is.null(del_group)) {
        del_group <- r
        del_start <- ref_pos
      } else {
        del_group <- c(del_group, r)
      }
    } else {
      # close a deletion group
      if (!is.null(del_group)) {
        del_end <- del_start + length(del_group) - 1
        if (length(del_group) == 1) {
          changes <- c(changes, paste0("p.", aa3[del_group], del_start, "del"))
        } else {
          changes <- c(changes, paste0("p.", aa3[del_group[1]], del_start, "_",
                                       aa3[del_group[length(del_group)]], del_end, "del"))
        }
        del_group <- NULL
        del_start <- NULL
      }
    }
    
    i <- i + 1
  }
  
  # any deletion group left open at alignment end
  if (!is.null(del_group)) {
    del_end <- del_start + length(del_group) - 1
    if (length(del_group) == 1) {
      changes <- c(changes, paste0("p.", aa3[del_group], del_start, "del"))
    } else {
      changes <- c(changes, paste0("p.", aa3[del_group[1]], del_start, "_",
                                   aa3[del_group[length(del_group)]], del_end, "del"))
    }
  }
  
  # Terminal insertions caught by alignment may have already been added.
  # Also add pure terminal insertion heuristics (in case pairwiseAlignment didn't create them)
  if (n_var > n_ref && startsWith(var_seq, ref_seq)) {
    ins_seq <- substring(var_seq, n_ref + 1)
    ins_3 <- paste(sapply(strsplit(ins_seq, "")[[1]], function(x) aa3[x]), collapse = "")
    last_res <- substr(ref_seq, n_ref, n_ref)
    changes <- c(changes, paste0("p.", aa3[last_res], n_ref, "_Ins", ins_3))
  }
  
  if (n_var > n_ref && endsWith(var_seq, ref_seq)) {
    ins_seq <- substring(var_seq, 1, n_var - n_ref)
    ins_3 <- paste(sapply(strsplit(ins_seq, "")[[1]], function(x) aa3[x]), collapse = "")
    first_res <- substr(ref_seq, 1, 1)
    # avoid duplicate when internal insertion between 1 and 2 already present
    second_res <- if (n_ref >= 2) substr(ref_seq, 2, 2) else NA_character_
    pattern_internal_12 <- NA
    if (!is.na(second_res)) {
      pattern_internal_12 <- paste0("^p\\.", aa3[first_res], "1_", aa3[second_res], "2ins")
    }
    if (is.na(pattern_internal_12) || !any(grepl(pattern_internal_12, changes, ignore.case = TRUE))) {
      changes <- c(changes, paste0("p.0_", aa3[first_res], "1ins", ins_3))
    }
  }
  
  # Terminal deletions (when variant shorter AND alignment leaves gaps at ends)
  # Use first/last aligned reference positions computed earlier:
  if (!is.na(first_ref_aligned_pos) && first_ref_aligned_pos > 1) {
    # N-terminal truncation was not explicitly requested, but if it occurs we report it
    del_start <- 1
    del_end <- first_ref_aligned_pos - 1
    changes <- c(changes, paste0("p.", aa3[substr(ref_seq,del_start,del_start)], del_start, "_",
                                 aa3[substr(ref_seq,del_end,del_end)], del_end, "del"))
  }
  if (!is.na(last_ref_aligned_pos) && last_ref_aligned_pos < n_ref) {
    del_start <- last_ref_aligned_pos + 1
    del_end <- n_ref
    changes <- c(changes, paste0("p.", aa3[substr(ref_seq,del_start,del_start)], del_start, "_",
                                 aa3[substr(ref_seq,del_end,del_end)], del_end, "del"))
  }
  
  paste(unique(changes), collapse = ";")
}

