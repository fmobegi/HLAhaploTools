# ~~~~~~~~~~~~~~~~~~~~
# HLAhaploTools()
# ~~~~~~~~~~~~~~~~~~~~
#' Run Full HLA Typing Pipeline with Optional Visualization
#'
#' Master wrapper that loads, cleans, optionally trims & decodes HLA typing data,
#' computes allele summaries, and runs EM‐based haplotype inference with plotting.
#'
#' @param filepath Path to the input HLA typing file (.csv, .tsv, .txt, .xlsx).
#' @param trim Logical. If TRUE, apply `trim_hla_results()` (default = FALSE).
#' @param resolution Integer. Trimming resolution (fields) for `trim_hla_results()` (default = 3).
#' @param isfamily Logical or NULL. If NULL (default), auto-detects data type. If TRUE, treats as family data.
#' @param mac Logical. If TRUE, decode multi-allele codes (default = FALSE).
#' @param plot_freq Logical. If TRUE, display allele frequency plot (default = FALSE).
#' @param plot_count Logical. If TRUE, display unique allele count plot (default = FALSE).
#' @param plot_diversity Logical. If TRUE, display diversity plot for specified gene(s) (default = FALSE).
#' @param plot_haplotypes Logical. If TRUE, display top haplotype frequency plot (default = FALSE).
#' @param gene Character vector. Gene(s) for diversity plotting (default = "A").
#' @param parallel Logical. If TRUE, use parallel processing when available (default = TRUE).
#' @param n_workers Integer. Number of workers for parallel processing (default = NULL, auto-detect).
#' @param sheet Excel sheet name or index to read (default = NULL).
#' @param quiet Logical. If TRUE, suppress status messages (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{typing_data}{Decoded, formatted HLA allele calls}
#'   \item{allele_frequencies}{Data frame of allele frequencies}
#'   \item{haplotype_frequencies}{Tibble of inferred haplotype frequencies}
#'   \item{posteriors}{List of full posterior diplotype distributions}
#'   \item{top_diplotypes}{Tibble of top diplotype call per subject}
#'   \item{loci_used}{Character vector of loci included in analysis}
#'   \item{convergence}{Logical indicating if EM algorithm converged}
#' }
#' @export
#'
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map
#' @importFrom tibble tibble
#' @importFrom dplyr mutate distinct slice_max relocate select
HLAhaploTools <- function(filepath,
                          trim = FALSE,
                          resolution = 3,
                          isfamily = NULL,
                          mac = TRUE,
                          plot_freq = FALSE,
                          plot_count = FALSE,
                          plot_diversity = FALSE,
                          plot_haplotypes = FALSE,
                          gene = "A",
                          parallel = TRUE,
                          n_workers = NULL,
                          sheet = NULL,
                          quiet = FALSE) {
  if (!quiet) {
    message("\n🧬 HLAhaploTools Analysis Pipeline")
  }

  # validate parameters and report which are set
  if (!quiet) {
    message("\n🔧 Parameters:")
    message(sprintf("\tTrim: %s (default resolution: %d)", trim, resolution))
    message(sprintf(
      "\tFamily data: %s",
      ifelse(is.null(isfamily), "Auto-detect", isfamily)
    ))
    message(sprintf("\tDecode MAC: %s", mac))
    message(sprintf("\tPlot allele frequency: %s", plot_freq))
    message(sprintf("\tPlot unique allele count: %s", plot_count))
    message(sprintf("\tPlot diversity: %s (gene: %s)", plot_diversity, gene))
    message(sprintf("\tPlot haplotypes: %s", plot_haplotypes))
    message(sprintf(
      "\tParallel processing: %s (workers: %s)",
      parallel,
      ifelse(is.null(n_workers), "Auto-detect", n_workers)
    ))
  }

  # Step 1: Load data file
  if (!quiet) {
    message("\n📁 Step 1: Loading HLA typing file...")
  }
  df_raw <- tryCatch(
    load_typing_data(filepath = filepath, sheet = sheet, quiet = quiet),
    error = function(e) {
      stop("\t❌ File load error: ", e$message)
    }
  )

  # Step 2: Detect or use provided data type info
  detect_result <- NULL
  if (is.null(isfamily)) {
    if (!quiet) message("\n🔍 Step 2: Auto-detecting data format...")
    detect_result <- detect_data_type(df_raw, quiet = quiet)
    family_data_val <- detect_result$is_family
  } else {
    family_data_val <- isfamily
    if (!quiet) {
      message(
        "\nℹ️ Data format explicitly provided as ",
        ifelse(isfamily, "'FAMILY STUDY'", "'REGULAR TYPING'"),
        "; skipping auto-detection."
      )
    }
    # Still run detection internally but quietly to get IDs and columns for downstream use
    detect_result <- detect_data_type(df_raw, quiet = TRUE)
  }

  # Assign relevant ID columns from detection
  family_col <- if (!is.null(detect_result$family_col)) {
    detect_result$family_col
  } else {
    NULL
  }
  member_col <- if (!is.null(detect_result$member_col)) {
    detect_result$member_col
  } else {
    NULL
  }
  id_cols <- if (!is.null(detect_result$id_cols)) {
    detect_result$id_cols
  } else {
    NULL
  }

  # Validate data structure based on detected type
  if (family_data_val) {
    if (!quiet) message("\n✓ Data format detected: FAMILY STUDY")
    valid <- validate_family_data(df_raw,
      stop_if_invalid = TRUE,
      verbose = !quiet
    )
  } else {
    if (!quiet) message("\n✓ Data format detected: REGULAR TYPING")
    valid <- validate_regular_data(df_raw,
      stop_if_invalid = TRUE,
      verbose = !quiet
    )
  }

  if (!valid) {
    stop("Validation failed: data does not meet expected format.")
  }
  # Step 3: Reformat typing data
  if (!quiet) {
    message("\n🧹 Step 3: Reformatting data...")
  }
  df_formatted <- reformat_typing_data(df_raw, isfamilydata = family_data_val, quiet = quiet)
  qc_report <- attr(df_formatted, "qc")

  # Step 4: Check deleted alleles
  if (!quiet) {
    message("\n⚠️ Step 4: Checking for deleted alleles...")
  }
  df_deleted_alleles <- check_deleted_alleles(df_formatted)

  # Step 5: If mac=TRUE, decode MAC strings
  if (mac) {
    if (!quiet) {
      message("\n🔄 Step 5: Decoding MAC strings...")
    }
    df_decoded <- decode_classical_mac(df_formatted, quiet = quiet)
  } else {
    df_decoded <- df_formatted
    if (!quiet) {
      message("\n➡️ Step 5: Skipping MAC decoding (mac = FALSE)")
    }
  }

  # Step 6: Optional trimming
  if (trim) {
    if (!quiet) {
      message("\n✂️ Step 6: Trimming allele names to specified resolution...")
    }
    df_decoded <- trim_hla_results(
      df_decoded,
      resolution = resolution,
      quiet = quiet
    )
  } else {
    if (!quiet) {
      message("\n➡️ Step 6: Skipping trimming (trim = FALSE)")
    }
  }

  # Step 7: Calculate allele frequencies
  if (!quiet) {
    message("\n📊 Step 7: Calculating allele frequencies...")
  }
  df_allele_freq <- calculate_hla_frequency(df_decoded, quiet = quiet)

  if (!quiet) {
    message("\t✅ Completed allele frequency summary.")

    # Get unique gene categories
    gene_names <- unique(df_allele_freq$gene)
    class1_genes <- intersect(c("A", "B", "C"), gene_names)
    class2_genes <- intersect(
      c(
        "DRB1",
        "DRB3",
        "DRB4",
        "DRB5",
        "DQA1",
        "DQB1",
        "DPA1",
        "DPB1"
      ),
      gene_names
    )
    nonclass_genes <- intersect(
      c("E", "F", "G", "H", "J", "K", "L", "HFE"),
      gene_names
    )
    mic_genes <- intersect(c("MICA", "MICB"), gene_names)

    message(sprintf(
      "\t📦 Genes: %d | Unique alleles: %d",
      length(gene_names),
      length(unique(df_allele_freq$allele))
    ))

    # Optional detailed gene reporting
    if (length(class1_genes) > 0) {
      message(sprintf(
        "\t✓ Class I genes: %s",
        paste(class1_genes, collapse = ", ")
      ))
    }
    if (length(class2_genes) > 0) {
      message(sprintf(
        "\t✓ Class II genes: %s",
        paste(class2_genes, collapse = ", ")
      ))
    }
    if (length(nonclass_genes) > 0) {
      message(sprintf(
        "\t✓ Non-classical genes: %s",
        paste(nonclass_genes, collapse = ", ")
      ))
    }
    if (length(mic_genes) > 0) {
      message(sprintf(
        "\t✓ MIC genes: %s",
        paste(mic_genes, collapse = ", ")
      ))
    }
  }

  # Step 8: Optional plots (allele freq, count, diversity) ...
  if (plot_freq) {
    if (!quiet) {
      message("\n🖼️  Step 7.1:  Generating allele frequency plot...")
    }
    print(plot_hla_allele_frequency(df_allele_freq, quiet = quiet))
  }

  if (plot_count) {
    if (!quiet) {
      message("\n📊  ️Step 7.2: Generating allele count plot...")
    }
    print(plot_hla_allele_count(df_allele_freq, quiet = quiet))
  }

  if (plot_diversity) {
    # Check if there are multiple populations to compare
    if ("population" %in% names(df_allele_freq) && length(unique(df_allele_freq$population)) > 1) {
      for (g in unique(gene)) {
        if (!quiet) {
          message(sprintf("\n🌐 Step 7.3: Generating diversity plot for gene: %s", g))
        }
        print(plot_hla_diversity(
          df_allele_freq,
          gene = g,
          quiet = quiet
        ))
      }
    } else {
      if (!quiet) {
        message("\n⚠️  Step 7.3: Diversity plot skipped: multiple populations required for comparison.")
      }
      message("\t   The current dataset contains only one population group.")
    }
  }

  # Step 9: Run haplotype inference with EM algorithm
  if (!quiet) {
    message("\n🔄 Step 9: Running haplotype inference with EM algorithm...")
  }

  haplotype_results <- infer_haplotypes(
    df = df_decoded,
    loci = NULL, # Auto-detect loci
    parallel = parallel,
    n_workers = n_workers,
    quiet = quiet,
    isfamily = family_data_val
  )

  hap_freqs <- haplotype_results$haplotype_frequencies
  posteriors <- haplotype_results$posteriors
  top_diplotypes <- haplotype_results$top_diplotypes

  if (!quiet) {
    message("\n✅ Haplotype inference complete!")
    message(
      sprintf(
        "\t📊 Found %d unique haplotypes (converged: %s, iterations: %d)",
        nrow(hap_freqs),
        ifelse(haplotype_results$convergence, "Yes", "No"),
        haplotype_results$iterations
      )
    )
  }

  # Step 10: Plot haplotype frequencies if requested
  if (plot_haplotypes) {
    if (!quiet) {
      message("\n📊 Generating haplotype frequency plot...")
    }
    print(plot_top_haplotypes(hap_freqs, quiet = quiet))
  }

  # Prepare result list
  result <- list(
    typing_data = df_decoded,
    allele_frequencies = df_allele_freq,
    haplotype_frequencies = hap_freqs,
    posteriors = posteriors,
    top_diplotypes = top_diplotypes,
    loci_used = haplotype_results$loci_used,
    convergence = haplotype_results$convergence
  )

  # Include deleted alleles if any rows present
  if (nrow(df_deleted_alleles) >= 1) {
    result$deleted_alleles <- df_deleted_alleles
  }

  if (!quiet) {
    message("\n🏁 HLAhaploTools analysis complete!")
  }

  return(result)
}
