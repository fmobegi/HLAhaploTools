# Example: clean column names for paired alleles
clean_hla_names <- function(df) {
  col_names <- colnames(df)

  corrected <- make.names(col_names, unique = TRUE)

  # Custom mapping: if Gene and Gene.1 exist, rename to Gene_1 and Gene_2
  for (gene in unique(sub("\\.\\d+$", "", corrected))) {
    matches <- grep(paste0("^", gene, "(\\.\\d+)?$"), corrected)
    if (length(matches) == 2) {
      corrected[matches] <- paste0(gene, c("_1", "_2"))
    }
  }

  colnames(df) <- corrected
  return(df)
}

hla_em_inference_dummy <- function(hla_df, loci = c("A", "B", "C", "DRB1", "DQA1", "DQB1"), quiet = FALSE) {
  if (!is.data.frame(hla_df) || nrow(hla_df) == 0) {
    stop("Input must be a non-empty data frame.")
  }

  # Step 0: filter out loci with excessive missingness
  usable_loci <- Filter(function(locus) {
    cols <- paste0(locus, c("_1", "_2"))
    all(cols %in% colnames(hla_df)) &&
      sum(rowSums(is.na(hla_df[, cols])) == 0) > 0
  }, loci)

  if (length(usable_loci) == 0) {
    stop("No loci with sufficient non-missing allele pairs available for EM inference.")
  }

  if (!quiet) {
    message(sprintf("\tℹ️ EM inference using loci: %s", paste(usable_loci, collapse = ", ")))
  }

  # Step 1: generate genotype pairs per locus
  genotype_pairs <- lapply(usable_loci, function(locus) {
    hla_df[, paste0(locus, c("_1", "_2"))]
  })
  names(genotype_pairs) <- usable_loci

  # Step 2: phased diplotypes per individual
  diplotype_sets <- lapply(seq_len(nrow(hla_df)), function(i) {
    allele_lists <- lapply(usable_loci, function(locus) {
      alleles <- na.omit(hla_df[i, paste0(locus, c("_1", "_2"))])
      if (length(alleles) == 2) {
        return(list(c(alleles[1], alleles[2])))
      } else {
        return(NULL)
      }
    })

    if (any(sapply(allele_lists, is.null))) {
      return(NULL) # skip individuals with incomplete typing
    }

    combos <- expand.grid(allele_lists)
    haplo_pairs <- list()
    for (j in seq_len(nrow(combos))) {
      h1 <- paste(mapply(`[`, combos[j, ], 1), collapse = "~")
      h2 <- paste(mapply(`[`, combos[j, ], 2), collapse = "~")
      haplo_pairs[[paste(h1, h2, sep = "|")]] <- 1
    }
    return(haplo_pairs)
  })

  diplotype_sets <- Filter(Negate(is.null), diplotype_sets)
  if (length(diplotype_sets) == 0) {
    stop("No complete diplotype sets could be constructed.")
  }

  # Step 3: extract unique haplotypes
  haplotypes <- unique(unlist(lapply(diplotype_sets, function(dset) {
    unlist(strsplit(names(dset), "\\|"))
  })))

  if (length(haplotypes) == 0) {
    stop("No valid haplotypes found.")
  }

  # Step 4: initialize frequencies
  haplo_freqs <- setNames(rep(1 / length(haplotypes), length(haplotypes)), haplotypes)

  # Step 5: EM algorithm
  max_iter <- 100
  epsilon <- 1e-5
  prev_ll <- -Inf

  for (iter in seq_len(max_iter)) {
    freq_counts <- setNames(rep(0, length(haplotypes)), haplotypes)
    ll <- 0

    for (diplo_set in diplotype_sets) {
      probs <- sapply(names(diplo_set), function(key) {
        haps <- strsplit(key, "\\|")[[1]]
        freq1 <- haplo_freqs[haps[1]]
        freq2 <- haplo_freqs[haps[2]]
        prob <- freq1 * freq2 * (1 + as.integer(haps[1] != haps[2]))
        if (is.na(prob)) 0 else prob
      })

      total <- sum(probs)
      if (!is.finite(total) || total == 0) next

      ll <- ll + log(total)
      norm_probs <- probs / total

      for (k in seq_along(norm_probs)) {
        pair <- strsplit(names(norm_probs)[k], "\\|")[[1]]
        p <- norm_probs[k]
        if (!is.na(p)) {
          freq_counts[pair[1]] <- freq_counts[pair[1]] + p
          freq_counts[pair[2]] <- freq_counts[pair[2]] + p
        }
      }
    }

    if (sum(freq_counts) == 0) {
      warning("No valid updates in EM step; stopping early.")
      break
    }

    freq_counts <- freq_counts / sum(freq_counts)
    delta <- max(abs(freq_counts - haplo_freqs), na.rm = TRUE)

    haplo_freqs <- freq_counts
    if (is.na(delta) || delta < epsilon) break
  }

  return(data.frame(Haplotype = names(haplo_freqs), Frequency = haplo_freqs, row.names = NULL))
}
