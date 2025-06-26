# ~~~~~~~~~~~~~~~~~~~~
# standardize_colnames()
# ~~~~~~~~~~~~~~~~~~~~
#' Standardize HLA Genotype Column Names
#'
#' Converts genotype column names (e.g. A.1) to A_1 format and applies consistent lowercase styling.
#'
#' @param df A raw HLA genotype data frame.
#' @return A cleaned data frame with standardized column names.
#' @export
#' @importFrom janitor clean_names
standardize_colnames <- function(df) {
  message("\t🔧 Standardizing genotype column names...")
  names(df) <- gsub("\\.", "_", names(df))
  janitor::clean_names(df, case = "none")
}

# ~~~~~~~~~~~~~~~~~~~~
# extract_loci()
# ~~~~~~~~~~~~~~~~~~~~
#' Extract HLA Locus Names from Genotype Table
#'
#' Identifies unique locus names based on *_1 and *_2 suffixes.
#'
#' @param df A data frame with standardized column names.
#' @return A character vector of distinct locus names (e.g. A, B, DRB1).
#' @export
extract_loci <- function(df) {
  message("\t🔬 Extracting HLA loci from standardized column names...")
  unique(gsub("_[12]$", "", names(df)[grepl("_[12]$", names(df))]))
}

# ~~~~~~~~~~~~~~~~~~~~
# collapse_genotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Collapse Allele Pairs Into Per-Locus Genotypes
#'
#' Combines paired allele fields (e.g. A_1 and A_2) into a string like A*01/A*02.
#'
#' @param df A data frame with per-locus allele columns.
#' @param loci A character vector of locus names to collapse.
#' @return A tibble with ID and one genotype column per gene.
#' @export
#' @importFrom dplyr bind_cols
#' @importFrom tibble tibble
collapse_genotypes <- function(df, loci) {
  message("\t🧬 Collapsing paired alleles into biallelic genotypes...")
  geno <- lapply(loci, function(gene) {
    a1 <- df[[paste0(gene, "_1")]]
    a2 <- df[[paste0(gene, "_2")]]
    val <- ifelse(is.na(a1) | is.na(a2),
      NA,
      paste(pmin(a1, a2), pmax(a1, a2), sep = "/")
    )
    tibble::tibble(!!gene := val)
  }) |> dplyr::bind_cols()
  dplyr::bind_cols(ID = df$ID, geno)
}

# ~~~~~~~~~~~~~~~~~~~~
# enumerate_diplotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' Exhaustive Diplotype Enumeration
#'
#' Returns all possible phased haplotype pairs consistent with a multilocus genotype row.
#'
#' @param geno_row A single row from collapsed genotype table.
#' @return A list of diplotypes with `hap1` and `hap2`.
#' @export
#' @importFrom utils combn
enumerate_diplotypes <- function(geno_row) {
  message("\t🔄 Generating exhaustive diplotype grid...")
  alleles_by_locus <- lapply(geno_row[-1], function(g) strsplit(g, "/", fixed = TRUE)[[1]])
  haplo_grid <- expand.grid(lapply(alleles_by_locus, combn, m = 1, simplify = FALSE)) |> as.data.frame()
  apply(haplo_grid, 1, function(h1) {
    h2 <- mapply(setdiff, alleles_by_locus, h1, SIMPLIFY = FALSE)
    list(
      hap1 = paste(unlist(h1), collapse = "~"),
      hap2 = paste(unlist(h2), collapse = "~")
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fast_diplotype_generator()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Fast Diplotype Generator (Non-Exhaustive)
#'
#' Infers exactly two phased diplotypes per individual by resolving alleles in input order.
#'
#' @param geno_row A single row of multilocus genotype.
#' @return A list of two diplotypes: hap1/hap2 and its reverse.
#' @export
fast_diplotype_generator <- function(geno_row) {
  message("\t🚀 Generating two-phase diplotypes (fast)...")
  loci <- names(geno_row)[-1]
  alleles <- lapply(geno_row[loci], function(g) strsplit(g, "/", fixed = TRUE)[[1]])
  valid_loci <- which(lengths(alleles) == 2 & !sapply(alleles, function(a) any(is.na(a))))
  if (length(valid_loci) < 1) {
    return(NULL)
  }
  alleles <- alleles[valid_loci]
  hap1 <- sapply(alleles, `[[`, 1)
  hap2 <- sapply(alleles, `[[`, 2)
  list(
    list(hap1 = paste(hap1, collapse = "~"), hap2 = paste(hap2, collapse = "~")),
    list(hap1 = paste(hap2, collapse = "~"), hap2 = paste(hap1, collapse = "~"))
  )
}

# ~~~~~~~~~~~~~~~~~~~~
# run_em_haplo()
# ~~~~~~~~~~~~~~~~~~~~
#' EM Algorithm for HLA Haplotype Frequency Estimation (Serial)
#'
#' Performs single‐threaded Expectation–Maximization to learn haplotype frequencies
#' from candidate diplotypes.
#'
#' @param all_diplotypes A list of diplotype lists (hap1/hap2) for each individual.
#' @param epsilon Numeric. Convergence threshold (default = 1e-5).
#' @param max_iter Integer. Maximum EM iterations (default = 100).
#'
#' @return A named numeric vector of estimated haplotype frequencies (sorted descending).
#' @export
run_em_haplo <- function(all_diplotypes, epsilon = 1e-5, max_iter = 100) {
  message("\t🔄 Initializing serial EM for haplotype frequency estimation...")
  hap_pool <- unique(unlist(lapply(all_diplotypes, function(dl) unlist(dl))))
  hap_freqs <- setNames(rep(1 / length(hap_pool), length(hap_pool)), hap_pool)

  for (iter in seq_len(max_iter)) {
    message(sprintf("\t⏳ EM iteration %d...", iter))

    # E-step
    exp_counts <- setNames(numeric(length(hap_freqs)), names(hap_freqs))
    for (dips in all_diplotypes) {
      if (is.null(dips)) next
      liks <- sapply(dips, function(dp) {
        h1 <- dp$hap1
        h2 <- dp$hap2
        if (h1 == h2) hap_freqs[h1]^2 else 2 * hap_freqs[h1] * hap_freqs[h2]
      })
      liks <- liks / sum(liks)
      for (j in seq_along(dips)) {
        h1 <- dips[[j]]$hap1
        h2 <- dips[[j]]$hap2
        exp_counts[h1] <- exp_counts[h1] + liks[j]
        exp_counts[h2] <- exp_counts[h2] + liks[j]
      }
    }

    # M-step & convergence
    new_freqs <- exp_counts / sum(exp_counts)
    diff <- max(abs(new_freqs - hap_freqs))
    message(sprintf("\t    Max freq change: %.2e", diff))
    if (diff < epsilon) {
      message(sprintf("\t✅ Converged after %d iterations.", iter))
      hap_freqs <- new_freqs
      break
    }
    hap_freqs <- new_freqs
  }

  message("\t📊 EM complete. Sorting haplotype frequencies.")
  sort(hap_freqs, decreasing = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run_em_haplo_multicore()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Parallel EM Algorithm for HLA Haplotype Frequencies
#'
#' Accelerates the EM E-step over individuals using future.apply parallelism.
#'
#' @inheritParams run_em_haplo
#' @export
#' @importFrom future.apply future_lapply
run_em_haplo_multicore <- function(all_diplotypes, epsilon = 1e-5, max_iter = 100) {
  message("\t🔄 Initializing parallel EM for haplotype frequency estimation...")
  hap_pool <- unique(unlist(lapply(all_diplotypes, function(dl) {
    if (is.null(dl)) {
      return(NULL)
    }
    unlist(lapply(dl, function(dp) c(dp$hap1, dp$hap2)))
  })))
  hap_freqs <- setNames(rep(1 / length(hap_pool), length(hap_pool)), hap_pool)

  for (iter in seq_len(max_iter)) {
    message(sprintf("\t⏳ Parallel EM iteration %d...", iter))

    # Parallel E-step
    exp_counts <- setNames(numeric(length(hap_freqs)), names(hap_freqs))
    partials <- future.apply::future_lapply(all_diplotypes, function(dips) {
      if (is.null(dips)) {
        return(NULL)
      }
      liks <- sapply(dips, function(dp) {
        h1 <- dp$hap1
        h2 <- dp$hap2
        if (h1 == h2) hap_freqs[h1]^2 else 2 * hap_freqs[h1] * hap_freqs[h2]
      })
      liks <- liks / sum(liks)
      local_counts <- setNames(numeric(length(hap_freqs)), names(hap_freqs))
      for (j in seq_along(dips)) {
        h1 <- dips[[j]]$hap1
        h2 <- dips[[j]]$hap2
        local_counts[h1] <- local_counts[h1] + liks[j]
        local_counts[h2] <- local_counts[h2] + liks[j]
      }
      local_counts
    })

    # Aggregate and M-step
    for (cnt in partials) {
      if (!is.null(cnt)) exp_counts <- exp_counts + cnt
    }
    new_freqs <- exp_counts / sum(exp_counts)
    diff <- max(abs(new_freqs - hap_freqs))
    message(sprintf("\t    Max freq change: %.2e", diff))

    if (diff < epsilon) {
      message(sprintf("\t✅ Converged after %d parallel iterations.", iter))
      hap_freqs <- new_freqs
      break
    }
    hap_freqs <- new_freqs
  }

  message("\t📊 Parallel EM complete. Sorting haplotype frequencies.")
  sort(hap_freqs, decreasing = TRUE)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign_HLA_posteriors()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Assign Posterior Diplotypes (Full &/or Top) Per Individual
#'
#' Computes posterior probabilities for all diplotypes per subject, deduplicates
#' reverse-equivalent pairs, and optionally returns both the full list and the
#' single top call.
#'
#' @param diplotype_list A list (length = n subjects) of diplotype candidates.
#'   Each element is a list of lists with fields `hap1` and `hap2`.
#' @param hap_freqs A named numeric vector of inferred haplotype frequencies.
#' @param ids Optional vector of subject IDs (length = n). Defaults to 1:n.
#' @param output Character. Which result to return:
#'   - `"both"`: return both full posterior list and top diplotypes (default)
#'   - `"full"`: return only the full posterior list
#'   - `"top"`: return only the top diplotypes tibble
#'
#' @return Depending on `output`:
#'   - A **list** with `posteriors` (list of tibbles) and `top_diplotypes` (tibble),
#'   - Or solely the `posteriors` list,
#'   - Or solely the `top_diplotypes` tibble.
#' @export
#'
#' @importFrom purrr map map_dfr
#' @importFrom tibble tibble
#' @importFrom dplyr mutate distinct slice_max relocate select desc
assign_HLA_posteriors <- function(diplotype_list,
                                  hap_freqs,
                                  ids = NULL,
                                  output = "both") {
  output <- match.arg(output, c("both", "full", "top"))
  n <- length(diplotype_list)
  if (is.null(ids)) ids <- seq_len(n)

  message("🧠 Computing posterior probabilities for each individual...")
  post_list <- purrr::map(seq_len(n), function(i) {
    dips <- diplotype_list[[i]]
    if (is.null(dips) || length(dips) == 0) {
      return(NULL)
    }

    # 1) Compute likelihoods and normalize
    probs <- sapply(dips, function(dp) {
      h1 <- dp$hap1
      h2 <- dp$hap2
      if (h1 == h2) hap_freqs[h1]^2 else 2 * hap_freqs[h1] * hap_freqs[h2]
    })
    total <- sum(probs)
    if (total == 0) {
      return(NULL)
    }
    post <- probs / total

    # 2) Build posterior tibble
    tibble::tibble(
      Haplotype1 = sapply(dips, `[[`, "hap1"),
      Haplotype2 = sapply(dips, `[[`, "hap2"),
      Posterior  = round(post, 6)
    ) |> dplyr::arrange(dplyr::desc(Posterior))
  })

  message("🎯 Extracting top diplotype per individual (deduplicated)...")
  top_tbl <- purrr::map_dfr(seq_len(n), function(i) {
    df <- post_list[[i]]
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    df %>%
      dplyr::mutate(
        h1  = pmin(Haplotype1, Haplotype2),
        h2  = pmax(Haplotype1, Haplotype2),
        Key = paste(h1, h2, sep = " || ")
      ) %>%
      dplyr::distinct(Key, .keep_all = TRUE) %>%
      dplyr::slice_max(order_by = Posterior, n = 1) %>%
      dplyr::mutate(ID = ids[i]) %>%
      dplyr::relocate(ID) %>%
      dplyr::select(-c(h1, h2, Key, Posterior))
  })

  # Return according to 'output'
  if (output == "full") {
    return(post_list)
  }
  if (output == "top") {
    return(top_tbl)
  }
  if (output == "both") {
    return(list(posteriors = post_list, top_diplotypes = top_tbl))
  }
}

# ~~~~~~~~~~~~~~~~~~~~
# choose_diplotype_expander()
# ~~~~~~~~~~~~~~~~~~~~
#' Select Diplotype Expansion Strategy Based on Core Count
#'
#' Chooses a fast or exhaustive diplotype generator function based on detected CPU cores.
#'
#' @return A function: either `fast_diplotype_generator` or `enumerate_diplotypes`.
#' @export
choose_diplotype_expander <- function() {
  n_cores <- parallel::detectCores(logical = TRUE)
  if (is.na(n_cores) || n_cores < 2) {
    message("\t💡 Not enough cores. Falling back to exhaustive diplotype enumerator...")
    return(enumerate_diplotypes)
  } else {
    message("\t🚀 Using fast diplotype generator (multi-core detected).")
    return(fast_diplotype_generator)
  }
}

# ~~~~~~~~~~~~~~~~~~~~
# choose_em_engine()
# ~~~~~~~~~~~~~~~~~~~~
#' Select EM Algorithm Engine Based on System Resources
#'
#' Chooses between single-threaded or parallel EM based on available CPU cores.
#'
#' @return A function: either `run_em_haplo` or `run_em_haplo_multicore`.
#' @export
choose_em_engine <- function() {
  n_cores <- parallel::detectCores(logical = TRUE)
  if (is.na(n_cores) || n_cores < 2) {
    message("\t💡 Not enough cores. Running EM in serial mode.")
    return(run_em_haplo)
  } else {
    message("\t🚀 Running parallel EM inference.")
    return(run_em_haplo_multicore)
  }
}

# ~~~~~~~~~~~~~~~~~~~~
# infer_haplotypes()
# ~~~~~~~~~~~~~~~~~~~~
#' EM-Based Inference of Multilocus HLA Haplotypes
#'
#' Runs EM algorithm on HLA genotype data to estimate haplotype frequencies
#' and assign most likely diplotype per individual.
#'
#' @param df A data frame of unphased HLA genotypes. Column names must follow <locus>_1 and <locus>_2 structure.
#' @param max_iter Maximum number of EM iterations. Default is 100.
#' @param epsilon Convergence threshold for EM algorithm. Default is 1e-5.
#'
#' @return A list with:
#' \describe{
#'   \item{frequencies}{A tibble of haplotypes and their estimated frequencies}
#'   \item{posteriors}{A list of posterior diplotype probabilities per individual}
#'   \item{top_diplotypes}{A tibble with one inferred diplotype per individual}
#' }
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate relocate select distinct slice_max

infer_haplotypes <- function(df, max_iter = 100, epsilon = 1e-5) {
  message("🧪 EM-Based Haplotype Inference Starting...")

  message("\t🧹 EM Step 1: Standardizing genotype format...")
  df_clean <- standardize_colnames(df)

  message("\t🧬 EM Step 2: Extracting loci and collapsing allele pairs...")
  loci <- extract_loci(df_clean)
  geno_df <- collapse_genotypes(df_clean, loci)

  message("\t⚙️ EM Step 3: Selecting diplotype generator and EM engine...")
  expand_fn <- choose_diplotype_expander()
  em_engine <- choose_em_engine()

  message("\t🔁 EM Step 4: Expanding diplotypes for each individual...")
  future::plan(future::multisession)
  diplotype_list <- furrr::future_map(
    1:nrow(geno_df),
    ~ expand_fn(geno_df[.x, ]),
    .progress = TRUE
  )

  message("\t📈 EM Step 5: Running EM to estimate haplotype frequencies...")
  hap_vector <- em_engine(diplotype_list, epsilon = epsilon, max_iter = max_iter)
  hap_freqs <- tibble::tibble(
    Haplotype = names(hap_vector),
    Frequency = round(as.numeric(hap_vector), 6)
  )

  message("\t🧠 EM Step 6: Computing posterior diplotype probabilities...")
  posterior_list <- compute_posteriors(diplotype_list, hap_vector)

  message("\t🎯 EM Step 7: Extracting top-scoring diplotype per individual...")
  top_diplos <- purrr::map_dfr(seq_along(posterior_list), function(i) {
    df_post <- posterior_list[[i]]
    if (is.null(df_post) || nrow(df_post) == 0) {
      return(NULL)
    }

    df_post %>%
      mutate(
        h1 = pmin(Haplotype1, Haplotype2),
        h2 = pmax(Haplotype1, Haplotype2),
        Key = paste(h1, h2, sep = " || ")
      ) %>%
      distinct(Key, .keep_all = TRUE) %>%
      slice_max(order_by = Posterior, n = 1) %>%
      mutate(ID = geno_df$ID[i]) %>%
      relocate(ID) %>%
      dplyr::select(-c(h1, h2, Key, Posterior))
  })

  message("\t✅ IE Inference complete.")
  return(list(
    frequencies = hap_freqs,
    posteriors = posterior_list,
    top_diplotypes = top_diplos
  ))
}
