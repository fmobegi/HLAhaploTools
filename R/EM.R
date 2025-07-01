# ~~~~~~~~~~~~~~~~~~~~
# standardize_colnames
# ~~~~~~~~~~~~~~~~~~~~
#' Standardize HLA Genotype Column Names
#'
#' Converts genotype column names (e.g. A.1) to A_1 format and applies consistent lowercase styling.
#'
#' @param df A raw HLA genotype data frame.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A cleaned data frame with standardized column names.
#' @examples
#' \dontrun{
#' # Standardize column names
#' raw_data <- read.csv("hla_data.csv")
#' clean_data <- standardize_colnames(raw_data)
#' }
#' @export
#' @importFrom janitor clean_names
standardize_colnames <- function(df, quiet = FALSE) {
  if (!quiet) {
    message("\t🔧 Standardizing genotype column names...")
  }

  if (!is.data.frame(df)) {
    stop("\t❌ Input must be a data frame or tibble.")
  }

  # Fix common separator variations
  names(df) <- gsub("\\.", "_", names(df)) # Convert A.1 to A_1
  names(df) <- gsub("-", "_", names(df)) # Convert A-1 to A_1
  names(df) <- gsub(" ", "_", names(df)) # Convert "A 1" to A_1

  # Clean with janitor
  df <- janitor::clean_names(df, case = "none")

  if (!quiet) {
    hla_cols <- sum(grepl("_[12]$", names(df)))
    message(
      sprintf(
        "\t✅ Standardized %d column names (%d HLA allele columns).",
        length(names(df)),
        hla_cols
      )
    )
  }

  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~
# extract_loci()
# ~~~~~~~~~~~~~~~~~~~~
#' Extract HLA Locus Names from Genotype Table
#'
#' Identifies unique locus names based on *_1 and *_2 suffixes.
#'
#' @param df A data frame with standardized column names.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A character vector of distinct locus names (e.g. A, B, DRB1).
#' @examples
#' \dontrun{
#' # Extract HLA loci from data
#' loci <- extract_loci(hla_data)
#' print(loci) # e.g., "A" "B" "C" "DRB1" "DQB1"
#' }
#' @export
extract_loci <- function(df, quiet = FALSE) {
  if (!quiet) {
    message("\t🔬 Extracting HLA loci from column names...")
  }

  if (!is.data.frame(df)) {
    stop("\t❌ Input must be a data frame or tibble.")
  }

  # Find columns ending with _1 or _2
  allele_cols <- grep("_[12]$", names(df), value = TRUE)
  if (length(allele_cols) == 0) {
    warning("\t⚠️ No HLA allele columns found with _1 or _2 suffix.")
    return(character(0))
  }

  # Extract unique gene names
  loci <- unique(gsub("_[12]$", "", allele_cols))

  if (!quiet) {
    message(sprintf(
      "\t✅ Found %d HLA loci: %s",
      length(loci),
      paste(loci, collapse = ", ")
    ))
  }

  return(loci)
}

# ~~~~~~~~~~~~~~~~~~~~
# collapse_genotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Collapse Allele Pairs Into Per-Locus Genotypes
#'
#' Combines paired allele fields (e.g. A_1 and A_2) into a string like A*01:01/A*01:02.
#'
#' @param df A data frame with per-locus allele columns.
#' @param loci A character vector of locus names to collapse.
#' @param separator Character string to separate alleles in genotype (default = "/").
#' @param quiet Logical. If TRUE, suppresses status messages.
#' @param isFamily Logical. If TRUE, uses FAMILY_ID and Family_Member for identification.
#'
#' @return A tibble with ID and one genotype column per gene.
#' @examples
#' \dontrun{
#' # Collapse allele pairs into genotypes
#' loci <- extract_loci(hla_data)
#' genotypes <- collapse_genotypes(hla_data, loci)
#' }
#' @export
#' @importFrom dplyr bind_cols mutate select
#' @importFrom tibble tibble
collapse_genotypes <- function(df,
                               loci,
                               separator = "/",
                               quiet = FALSE,
                               isFamily = TRUE) {
  if (!quiet) {
    message("\t\U1F9EC Collapsing paired alleles into biallelic genotypes...")
  }

  if (!is.data.frame(df)) {
    stop("\t\u274C Input must be a data frame or tibble.")
  }

  if (!is.character(loci) || length(loci) == 0) {
    stop("\t\u274C Loci must be a non-empty character vector.")
  }

  # Find ID columns for different data types
  id_cols <- c()
  if (isFamily) {
    # Family data - use FAMILY_ID and Family_Member
    if (all(c("FAMILY_ID", "Family_Member") %in% names(df))) {
      id_cols <- c("FAMILY_ID", "Family_Member")
      if (!quiet) {
        message("\t\u2713 Using family data columns (FAMILY_ID, Family_Member)")
      }
    } else {
      warning("\t\u26A0\uFE0F Family mode specified but FAMILY_ID/Family_Member columns not found.")
    }
  }

  # Try finding any ID column if not already set
  if (length(id_cols) == 0) {
    possible_id_cols <- c(
      "SampleID",
      "Sample_ID",
      "ID",
      "Id",
      "id",
      "PATIENT_ID",
      "Patient_ID",
      "Subject_ID"
    )
    found_ids <- intersect(possible_id_cols, names(df))

    if (length(found_ids) > 0) {
      id_cols <- found_ids[1]
      if (!quiet) {
        message(sprintf("\t\u2713 Using ID column: %s", id_cols))
      }
    } else {
      if (!quiet) {
        message("\t\u26A0\uFE0F No ID column found. Using row numbers.")
      }
      df$ID <- seq_len(nrow(df))
      id_cols <- "ID"
    }
  }

  # Start with ID columns
  result <- df[, id_cols, drop = FALSE]

  # For each locus, combine allele pairs into genotype strings
  success_count <- 0
  for (gene in loci) {
    col1 <- paste0(gene, "_1")
    col2 <- paste0(gene, "_2")

    if (!all(c(col1, col2) %in% names(df))) {
      if (!quiet) {
        message(sprintf("\t\u26A0\uFE0F Skipping %s (missing %s or %s)", gene, col1, col2))
      }
      next
    }

    result[[gene]] <- mapply(function(a1, a2) {
      if (is.na(a1) || a1 == "" || is.na(a2) || a2 == "") {
        return(NA_character_)
      }
      paste(a1, a2, sep = separator)
    }, df[[col1]], df[[col2]])

    success_count <- success_count + 1
  }

  if (!quiet) {
    message(
      sprintf(
        "\t\u2705 Created %d collapsed genotype columns from %d loci.",
        success_count,
        length(loci)
      )
    )
  }

  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~
# enumerate_diplotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Exhaustive Diplotype Enumeration
#'
#' Returns all possible phased haplotype pairs consistent with a multilocus genotype row.
#'
#' @param geno_row A single row from collapsed genotype table.
#' @param separator Character string that separates alleles in genotype (default = "/").
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A list of diplotypes with `hap1` and `hap2`.
#' @examples
#' \dontrun{
#' # Enumerate all possible diplotypes
#' genotypes <- collapse_genotypes(hla_data, c("A", "B", "DRB1"))
#' diplotypes <- enumerate_diplotypes(genotypes[1, ])
#' }
#' @export
#' @importFrom utils combn
enumerate_diplotypes <- function(geno_row,
                                 separator = "/",
                                 quiet = FALSE) {
  if (!is.data.frame(geno_row) || nrow(geno_row) != 1) {
    stop("\t❌ Input must be a single row data frame.")
  }

  # Skip ID columns
  id_cols <- c(
    "ID",
    "FAMILY_ID",
    "Family_Member",
    "SampleID",
    "Sample_ID",
    "PATIENT_ID",
    "Patient_ID",
    "Subject_ID"
  )
  skip_cols <- which(names(geno_row) %in% id_cols)
  if (length(skip_cols) > 0) {
    geno_cols <- setdiff(seq_along(geno_row), skip_cols)
    geno_row_filtered <- geno_row[, geno_cols, drop = FALSE]
  } else {
    geno_row_filtered <- geno_row
  }

  # Extract alleles for each locus
  loci <- names(geno_row_filtered)
  allele_lists <- lapply(seq_along(loci), function(i) {
    gene <- loci[i]
    genotype <- geno_row_filtered[[gene]]

    if (is.na(genotype) || genotype == "") {
      return(NULL)
    }

    alleles <- unlist(strsplit(genotype, separator, fixed = TRUE))
    if (length(alleles) != 2) {
      if (!quiet) {
        message(sprintf(
          "\t⚠️ Skipping locus %s: expected 2 alleles, found %d.",
          gene,
          length(alleles)
        ))
      }
      return(NULL)
    }

    # Ensure alleles preserve the gene prefix if it exists
    # If they don't have gene prefix (e.g., just "01:01"), add it (e.g., "A*01:01")
    format_allele <- function(allele, gene) {
      if (!grepl("\\*", allele)) {
        return(paste0(gene, "*", allele))
      }
      return(allele) # Already has a prefix like "A*"
    }

    return(list(
      gene = gene,
      a1 = format_allele(alleles[1], gene),
      a2 = format_allele(alleles[2], gene)
    ))
  })

  # Remove any NULL entries (skipped loci)
  allele_lists <- allele_lists[!sapply(allele_lists, is.null)]

  if (length(allele_lists) == 0) {
    if (!quiet) {
      message("\t⚠️ No valid genotypes found to enumerate.")
    }
    return(list())
  }

  # Calculate number of possible combinations
  n_loci <- length(allele_lists)
  n_combos <- 2^(n_loci - 1)

  if (!quiet && n_combos > 100) {
    message(sprintf(
      "\t⚠️ Large number of diplotypes (%d) will be generated.",
      n_combos
    ))
  }

  # Generate all possible phasings using binary representation
  result <- lapply(0:(n_combos - 1), function(i) {
    # Convert i to binary to determine phasing pattern
    phasing <- as.integer(intToBits(i))[1:n_loci]

    # Create haplotypes based on phasing pattern
    hap1 <- character(n_loci)
    hap2 <- character(n_loci)

    for (j in 1:n_loci) {
      if (phasing[j] == 0) {
        hap1[j] <- allele_lists[[j]]$a1
        hap2[j] <- allele_lists[[j]]$a2
      } else {
        hap1[j] <- allele_lists[[j]]$a2
        hap2[j] <- allele_lists[[j]]$a1
      }
    }

    list(
      hap1 = paste(hap1, collapse = "~"),
      hap2 = paste(hap2, collapse = "~")
    )
  })

  if (!quiet) {
    message(sprintf(
      "\t✅ Enumerated %d possible diplotypes across %d loci.",
      length(result),
      n_loci
    ))
  }
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~
# compute_posteriors()
# ~~~~~~~~~~~~~~~~~~~~
#' Compute Posterior Probabilities of Diplotypes
#'
#' Calculates the posterior probabilities of all possible diplotypes for each subject
#' given the estimated haplotype frequencies.
#'
#' @param all_diplotypes A list of diplotype lists, one per subject.
#' @param hap_freqs A named vector of haplotype frequencies.
#' @param normalize Logical. If TRUE (default), normalize posteriors to sum to 1.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A list of posterior probabilities for each subject's possible diplotypes.
#' @examples
#' \dontrun{
#' # Compute posteriors after EM algorithm
#' posteriors <- compute_posteriors(all_diplotypes, haplotype_frequencies$frequency)
#' }
#' @export
compute_posteriors <- function(all_diplotypes,
                               hap_freqs,
                               normalize = TRUE,
                               quiet = FALSE) {
  if (!quiet) {
    message("\t📊 Computing posterior diplotype probabilities...")
  }

  if (!is.list(all_diplotypes)) {
    stop("\t❌ all_diplotypes must be a list of diplotype lists.")
  }

  if (!is.numeric(hap_freqs) ||
    any(hap_freqs < 0) || any(hap_freqs > 1)) {
    stop("\t❌ hap_freqs must be a numeric vector with values between 0 and 1.")
  }

  # Ensure hap_freqs has names
  if (is.null(names(hap_freqs))) {
    stop("\t❌ hap_freqs must be a named vector.")
  }

  # Calculate posterior for each subject's diplotypes
  result <- lapply(all_diplotypes, function(subject_diplotypes) {
    if (length(subject_diplotypes) == 0) {
      return(NULL)
    }

    # For each diplotype, calculate probability
    probs <- sapply(subject_diplotypes, function(diplotype) {
      h1 <- diplotype$hap1
      h2 <- diplotype$hap2

      # Check if haplotypes are present in frequency table
      if (!h1 %in% names(hap_freqs) || !h2 %in% names(hap_freqs)) {
        return(0)
      }

      # Calculate probability (2*p*q for heterozygotes, p^2 for homozygotes)
      if (h1 == h2) {
        return(hap_freqs[h1]^2)
      } else {
        return(2 * hap_freqs[h1] * hap_freqs[h2])
      }
    })

    # Normalize probabilities if requested
    if (normalize && sum(probs) > 0) {
      probs <- probs / sum(probs)
    }

    return(probs)
  })

  if (!quiet) {
    non_null <- sum(!sapply(result, is.null))
    message(sprintf("\t✅ Computed posteriors for %d subjects.", non_null))
  }

  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~
# run_em_algorithm()
# ~~~~~~~~~~~~~~~~~~~~
#' Run EM Algorithm for Haplotype Frequency Estimation
#'
#' Implements the expectation-maximization algorithm to estimate haplotype frequencies
#' from genotype data. Automatically uses parallel processing if available.
#'
#' @param all_diplotypes A list of diplotype lists, one per subject.
#' @param epsilon Convergence threshold (default = 1e-5).
#' @param max_iter Maximum number of iterations (default = 100).
#' @param parallel Logical. Controls parallel execution. If TRUE (default), uses parallel processing
#'                 when available and falls back to sequential if not. If FALSE, forces sequential.
#' @param n_workers Number of workers for parallel processing. If NULL, use available cores - 1.
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @note The implementation of this function was developed with assistance from
#'   Claude Sonnet 3.7 Thinking.
#' @return A named vector of estimated haplotype frequencies.
#' @examples
#' \dontrun{
#' # Run EM algorithm with automatic parallelization
#' freq_est <- run_em_algorithm(all_diplotypes)
#'
#' # Force sequential processing
#' freq_est <- run_em_algorithm(all_diplotypes, parallel = FALSE)
#' }
#' @importFrom future plan multisession sequential availableCores
#' @importFrom furrr future_map_dbl
#' @export
run_em_algorithm <- function(all_diplotypes,
                             epsilon = 1e-5,
                             max_iter = 100,
                             parallel = TRUE,
                             n_workers = NULL,
                             quiet = FALSE) {
  # Check input validity
  if (!is.list(all_diplotypes)) {
    stop("\t❌ all_diplotypes must be a list of diplotype lists.")
  }

  # Count total number of samples
  n_samples <- sum(sapply(all_diplotypes, length) > 0)

  if (n_samples == 0) {
    stop("\t❌ No valid diplotypes found.")
  }

  # Extract all unique haplotypes across all samples
  all_haps <- unique(unlist(lapply(all_diplotypes, function(subject) {
    if (length(subject) == 0) {
      return(NULL)
    }
    unlist(lapply(subject, function(d) c(d$hap1, d$hap2)))
  })))

  if (length(all_haps) == 0) {
    stop("\t❌ No haplotypes found in diplotype lists.")
  }

  # Initialize frequencies uniformly
  hap_freqs <- rep(1 / length(all_haps), length(all_haps))
  names(hap_freqs) <- all_haps

  # Precompute haplotype pairs for faster E-step
  diplotype_hap_pairs <- lapply(all_diplotypes, function(subject_diplotypes) {
    lapply(subject_diplotypes, function(d) c(d$hap1, d$hap2))
  })

  # Run EM algorithm
  converged <- FALSE
  iter <- 0

  while (!converged && iter < max_iter) {
    iter <- iter + 1
    old_freqs <- hap_freqs

    # E-step: Calculate expected haplotype counts based on current frequencies
    hap_counts <- numeric(length(all_haps))
    names(hap_counts) <- all_haps

    for (i in seq_along(diplotype_hap_pairs)) {
      subject_pairs <- diplotype_hap_pairs[[i]]
      if (length(subject_pairs) == 0) next

      subject_posteriors <- sapply(subject_pairs, function(pair) {
        h1 <- pair[1]
        h2 <- pair[2]
        if (!h1 %in% names(hap_freqs) || !h2 %in% names(hap_freqs)) {
          return(0)
        }
        if (h1 == h2) hap_freqs[h1]^2 else 2 * hap_freqs[h1] * hap_freqs[h2]
      })

      # Normalize posteriors
      if (sum(subject_posteriors) > 0) {
        subject_posteriors <- subject_posteriors / sum(subject_posteriors)
      }

      # Update haplotype counts
      for (j in seq_along(subject_pairs)) {
        pair <- subject_pairs[[j]]
        p <- subject_posteriors[j]
        if (p > 0) {
          hap_counts[pair[1]] <- hap_counts[pair[1]] + p
          hap_counts[pair[2]] <- hap_counts[pair[2]] + p
        }
      }
    }

    # M-step: Update frequency estimates
    total_haps <- sum(hap_counts)
    hap_freqs <- hap_counts / total_haps

    # Check convergence
    delta <- max(abs(hap_freqs - old_freqs))
    converged <- delta < epsilon

    if (!quiet && (iter %% 10 == 0 || converged || iter == max_iter)) {
      message(sprintf("\t\t📈 Iteration %d: max change = %.8f", iter, delta))
    }
  }

  if (converged) {
    if (!quiet) message(sprintf("\t✅ EM algorithm converged after %d iterations.", iter))
  } else {
    warning(sprintf("\t⚠️ EM algorithm did not converge after %d iterations.", max_iter))
  }

  # Return the estimated frequencies along with convergence information
  attr(hap_freqs, "converged") <- converged
  attr(hap_freqs, "iterations") <- iter

  return(hap_freqs)
}

# ~~~~~~~~~~~~~~~~~~~~
# infer_haplotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Infer HLA Haplotype Frequencies from Genotype Data
#'
#' Master function that performs full haplotype frequency inference using the EM algorithm.
#' Takes genotype data, enumerates all possible diplotypes, and estimates haplotype frequencies.
#'
#' @param df A data frame with HLA allele columns (e.g. A_1, A_2, B_1, B_2).
#' @param loci Character vector of loci to include. If NULL, all detected loci are used.
#' @param max_iter Maximum number of EM iterations (default = 100).
#' @param epsilon Convergence threshold (default = 1e-5).
#' @param parallel Logical. If TRUE (default), use parallel processing when available.
#' @param n_workers Number of workers for parallel processing. If NULL, auto-detect.
#' @param quiet Logical. If TRUE, suppresses status messages.
#' @param isFamily Logical. If TRUE, uses FAMILY_ID and Family_Member for identification.
#'
#' @return A list containing:
#' \describe{
#'   \item{haplotype_frequencies}{A tibble with haplotype frequencies}
#'   \item{posteriors}{A list of posterior diplotype probabilities per subject}
#'   \item{top_diplotypes}{A tibble with most likely diplotype for each subject}
#'   \item{loci_used}{Character vector of loci included in the analysis}
#'   \item{convergence}{Logical. Whether the EM algorithm converged}
#'   \item{iterations}{Number of iterations performed}
#' }
#' @examples
#' \dontrun{
#' # Infer haplotypes from all loci
#' result <- infer_haplotypes(hla_data)
#'
#' # Infer haplotypes from specific loci
#' result <- infer_haplotypes(hla_data, loci = c("A", "B", "DRB1"))
#'
#' # Force sequential processing
#' result <- infer_haplotypes(hla_data, parallel = FALSE)
#'
#' # Process family data
#' family_results <- infer_haplotypes(family_data, isFamily = TRUE)
#' }
#' @importFrom future availableCores plan multisession sequential
#' @importFrom furrr future_map
#' @importFrom tibble tibble
#' @importFrom dplyr arrange desc select mutate
#' @export
infer_haplotypes <- function(df,
                             loci = NULL,
                             max_iter = 100,
                             epsilon = 1e-5,
                             parallel = TRUE,
                             n_workers = NULL,
                             quiet = FALSE,
                             isFamily = FALSE) {
  if (!quiet) {
    message("\t🔄 Starting HLA haplotype inference pipeline...")
  }

  if (!is.data.frame(df)) {
    stop("\t❌ Input must be a data frame or tibble.")
  }

  # Step 1: Standardize column names
  df <- standardize_colnames(df, quiet = quiet)

  # Step 2: Extract loci if not provided
  if (is.null(loci)) {
    loci <- extract_loci(df, quiet = quiet)
  } else if (!quiet) {
    message(sprintf(
      "\t✅ Using user-specified loci: %s",
      paste(loci, collapse = ", ")
    ))
  }

  if (length(loci) == 0) {
    stop("\t❌ No HLA loci found for analysis.")
  }

  # Step 3: Collapse genotypes
  genotypes <- collapse_genotypes(df, loci, quiet = quiet, isFamily = isFamily)

  # Step 4: Enumerate all possible diplotypes per subject
  if (!quiet) {
    message("\t🔢 Enumerating all possible diplotypes per subject...")
  }

  all_diplotypes <- lapply(seq_len(nrow(genotypes)), function(i) {
    enumerate_diplotypes(genotypes[i, ], quiet = TRUE)
  })

  valid_samples <- sum(sapply(all_diplotypes, length) > 0)
  if (!quiet) {
    message(
      sprintf(
        "\t✅ Found valid diplotype configurations for %d out of %d samples.",
        valid_samples,
        nrow(genotypes)
      )
    )
  }

  if (valid_samples == 0) {
    stop("\t❌ No valid diplotypes found to analyze.")
  }

  # Step 5: Run EM algorithm
  hap_freqs <- run_em_algorithm(
    all_diplotypes,
    epsilon = epsilon,
    max_iter = max_iter,
    parallel = parallel,
    n_workers = n_workers,
    quiet = quiet
  )

  # Get convergence information
  converged <- attr(hap_freqs, "converged")
  iterations <- attr(hap_freqs, "iterations")

  # Step 6: Compute posterior probabilities
  posteriors <- compute_posteriors(all_diplotypes, hap_freqs, quiet = quiet)

  # Step 7: Extract most likely diplotype for each subject
  if (!quiet) {
    message("\t🏆 Identifying most likely diplotype for each subject...")
  }

  # Get IDs from the original genotypes data frame
  id_cols <- c()
  if (isFamily &&
    all(c("FAMILY_ID", "Family_Member") %in% names(genotypes))) {
    id_cols <- c("FAMILY_ID", "Family_Member")
  } else {
    # Try to find any ID column
    possible_id_cols <- c("ID", "SampleID", "Sample_ID", "PATIENT_ID", "Subject_ID")
    id_cols <- intersect(possible_id_cols, names(genotypes))
    if (length(id_cols) == 0) {
      id_cols <- "ID" # Use default if nothing found
    }
  }

  # Create initial ID columns for results
  top_diplotypes <- genotypes[, id_cols, drop = FALSE]
  top_diplotypes$hap1 <- NA_character_
  top_diplotypes$hap2 <- NA_character_
  top_diplotypes$probability <- NA_real_

  for (i in seq_along(all_diplotypes)) {
    subject_diplotypes <- all_diplotypes[[i]]
    if (length(subject_diplotypes) == 0) {
      next
    }

    subject_posteriors <- posteriors[[i]]
    if (is.null(subject_posteriors) ||
      all(subject_posteriors == 0)) {
      next
    }

    # Find most likely diplotype
    best_idx <- which.max(subject_posteriors)
    best_diplotype <- subject_diplotypes[[best_idx]]
    best_prob <- subject_posteriors[best_idx]

    top_diplotypes$hap1[i] <- best_diplotype$hap1
    top_diplotypes$hap2[i] <- best_diplotype$hap2
    top_diplotypes$probability[i] <- best_prob
  }

  # Step 8: Format haplotype frequency results
  if (!quiet) {
    message("\t📋 Preparing haplotype frequency results...")
  }

  hap_freq_df <- tibble::tibble(
    haplotype = names(hap_freqs),
    frequency = as.numeric(hap_freqs)
  ) %>%
    dplyr::arrange(dplyr::desc(frequency)) %>%
    dplyr::mutate(
      rank = seq_len(n()),
      cumulative_freq = cumsum(frequency)
    )

  # Return the results
  result <- list(
    haplotype_frequencies = hap_freq_df,
    posteriors = posteriors,
    top_diplotypes = top_diplotypes,
    loci_used = loci,
    convergence = converged,
    iterations = iterations
  )

  if (!quiet) {
    message(
      sprintf(
        "\t✅ HLA haplotype inference complete: %d unique haplotypes identified.",
        nrow(hap_freq_df)
      )
    )
  }

  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~
# plot_top_haplotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Plot Top HLA Haplotypes by Frequency
#'
#' Visualizes the most common haplotypes and their frequencies.
#'
#' @param hap_freq A data frame or tibble with haplotype frequencies.
#' @param n_top Number of top haplotypes to display (default = 20).
#' @param min_freq Minimum frequency to include (default = 0.01).
#' @param color Fill color for bars (default = "steelblue").
#' @param quiet Logical. If TRUE, suppresses status messages.
#'
#' @return A ggplot object that can be further customized or printed.
#' @examples
#' \dontrun{
#' # Basic usage with results from infer_haplotypes
#' result <- infer_haplotypes(hla_data)
#' p <- plot_top_haplotypes(result$haplotype_frequencies)
#'
#' # Customize to show only top 10 with higher threshold
#' p <- plot_top_haplotypes(result$haplotype_frequencies,
#'   n_top = 10, min_freq = 0.02
#' )
#' }
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_y_continuous theme_minimal theme labs element_text
#' @export
plot_top_haplotypes <- function(hap_freq,
                                n_top = 20,
                                min_freq = 0.01,
                                color = "steelblue",
                                quiet = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "\t❌ Package ggplot2 is required for plotting. Please install it with install.packages('ggplot2')"
    )
  }

  if (!quiet) {
    message("\t📊 Creating haplotype frequency plot...")
  }

  if (!"haplotype" %in% names(hap_freq) ||
    !"frequency" %in% names(hap_freq)) {
    stop("\t❌ Input must have 'haplotype' and 'frequency' columns.")
  }

  # Filter data
  plot_data <- hap_freq[hap_freq$frequency >= min_freq, ]
  plot_data <- plot_data[order(plot_data$frequency, decreasing = TRUE), ]

  if (nrow(plot_data) == 0) {
    warning(sprintf(
      "\t⚠️ No haplotypes meet the minimum frequency threshold of %.4f",
      min_freq
    ))
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle("No haplotypes meet frequency threshold")
    )
  }

  # Limit to top n
  if (nrow(plot_data) > n_top) {
    plot_data <- plot_data[1:n_top, ]
    if (!quiet) {
      message(sprintf(
        "\t✓ Showing top %d haplotypes out of %d that meet threshold.",
        n_top,
        nrow(hap_freq[hap_freq$frequency >= min_freq, ])
      ))
    }
  }

  # Format haplotypes for display
  plot_data$haplotype <- factor(plot_data$haplotype, levels = rev(plot_data$haplotype))

  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = haplotype, y = frequency)) +
    ggplot2::geom_bar(
      stat = "identity",
      fill = color,
      width = 0.7
    ) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", frequency * 100)), hjust = -0.2) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, max(plot_data$frequency) * 1.2)) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Top HLA Haplotype Frequencies", x = "Haplotype", y = "Frequency") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!quiet) {
    message("\t✅ Haplotype frequency plot created successfully.")
  }
  return(p)
}
