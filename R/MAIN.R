# ~~~~~~~~~~~~~~~~~~~~
# HLAhaploTools()
# ~~~~~~~~~~~~~~~~~~~~
#' Run Full HLA Typing Pipeline with Optional Visualization
#'
#' A master wrapper that loads HLA typing data from file, detects data format,
#' applies optional formatting and trimming, decodes multi-allele strings (MAC),
#' calculates allele-level summaries, performs haplotype inference via EM algorithm,
#' and optionally generates plots. Supports both family and registry/population-based studies.
#'
#' @param filepath Character. Path to the input HLA typing file (.csv, .tsv, .txt, .xlsx).
#' @param trim Logical. If TRUE, applies `trim_hla_results()` to standardize allele
#'   resolution. Default is FALSE.
#' @param resolution Integer. Desired resolution for allele trimming. Default is 3.
#' @param isfamily Logical or NULL. If NULL (default), auto-detects data format.
#'   If TRUE, treats input as family data.
#' @param mac Logical. If TRUE, decode multi-allele codes using `decode_classical_mac()`.
#'   Default is TRUE.
#' @param plot_freq Logical. If TRUE, generates allele frequency plot. Default is FALSE.
#' @param plot_count Logical. If TRUE, plots unique allele counts across loci.
#'   Default is FALSE.
#' @param plot_diversity Logical. If TRUE, plots diversity across populations for
#'   selected gene(s). Default is FALSE.
#' @param plot_haplotypes Logical. If TRUE, plots top inferred haplotype frequencies.
#'   Default is FALSE.
#' @param gene Character vector. Genes to use for diversity plotting. Default is "A".
#' @param parallel Logical. If TRUE, uses parallelization for EM haplotype inference.
#'   Default is TRUE.
#' @param n_workers Integer or NULL. Number of parallel workers. If NULL, uses
#'   `availableCores()`. Default is NULL.
#' @param sheet Character or integer. Sheet name or index to use when reading Excel
#'   files. Optional.
#' @param quiet Logical. If TRUE, suppresses progress and status messaging.
#'   Default is FALSE.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{typing_data},{Cleaned and optionally trimmed HLA typing data}
#'   \item{allele_frequencies},{Data frame of allele frequency summaries}
#'   \item{haplotype_frequencies},{Tibble containing inferred haplotype distributions}
#'   \item{posteriors},{List of posterior diplotype probabilities per subject}
#'   \item{top_diplotypes},{Tibble summarizing most likely diplotype assignment per subject}
#'   \item{loci_used},{Vector of HLA loci used in haplotype inference}
#'   \item{convergence},{Logical indicating convergence status of EM algorithm}
#'   \item{deleted_alleles},{Optional. Tibble of alleles flagged as deleted in IMGT/HLA reference}
#'   \item{family_segregation},{Optional. Output of segregation analysis if input is family-based}
#' }
#'
#' @details
#' This function provides a comprehensive pipeline for HLA typing analysis including:
#' \itemize{
#'   \item Data loading and format detection
#'   \item Multi-allele code (MAC) decoding
#'   \item Allele frequency calculations
#'   \item Haplotype inference using EM algorithm
#'   \item Family segregation analysis (for family data)
#'   \item Optional visualization plots
#' }
#'
#' The function automatically detects whether the input data represents family studies
#' or population/registry data and applies appropriate analysis methods.
#'
#' @examples
#' \dontrun{
#' # Basic usage with family data
#' results <- HLAhaploTools(
#'    filepath = "inst/extdata/family_typing_data.txt",
#'    isfamily = TRUE,
#'    plot_freq = TRUE,
#'    plot_haplotypes = TRUE
#' )
#'
#' # Advanced usage with all options
#' results <- HLAhaploTools(
#'    filepath = "inst/extdata/family_typing_data.txt",
#'    trim = FALSE,
#'    resolution = 3,
#'    isfamily = TRUE,
#'    plot_freq = TRUE,
#'    plot_count = TRUE,
#'    plot_diversity = TRUE,
#'    plot_haplotypes = TRUE,
#'    gene = "A",
#'    parallel = TRUE,
#'    n_workers = 8,
#'    quiet = FALSE,
#'    mac = TRUE
#' )
#'
#' # Access results
#' typing_data <- results$typing_data
#' allele_freq <- results$allele_frequencies
#' haplotypes <- results$haplotype_frequencies
#' }
#'
#' @seealso
#' \code{\link{load_typing_data}}, \code{\link{reformat_typing_data}},
#' \code{\link{decode_classical_mac}}, \code{\link{infer_haplotypes}}
#'
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map
#' @importFrom tibble tibble
#' @importFrom dplyr mutate distinct slice_max relocate select
#'
#' @export
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
   if (!quiet) message("\n🧬 HLAhaploTools Analysis Pipeline")

   # Parameters
   if (!quiet) {
      message("\n🔧 Parameters:")
      message(sprintf("   Trim: %s (default resolution: %d)", trim, resolution))
      message(sprintf("   Family data: %s", ifelse(is.null(isfamily), "Auto-detect", isfamily)))
      message(sprintf("   Decode MAC: %s", mac))
      message(sprintf("   Plot allele frequency: %s", plot_freq))
      message(sprintf("   Plot unique allele count: %s", plot_count))
      message(sprintf("   Plot diversity: %s (gene: %s)", plot_diversity, gene))
      message(sprintf("   Plot haplotypes: %s", plot_haplotypes))
      message(sprintf(
         "   Parallel processing: %s (workers: %s)",
         parallel,
         ifelse(is.null(n_workers), "Auto-detect", n_workers)
      ))
   }

   # Step 1: Load data file
   if (!quiet) message("\n📁 Step 1: Loading HLA typing file...")
   df_raw <- tryCatch(
      load_typing_data(filepath = filepath, sheet = sheet, quiet = quiet),
      error = function(e) stop("   ❌ File load error: ", e$message)
   )

   # Step 2: Detect or use provided data type info
   if (is.null(isfamily)) {
      if (!quiet) message("\n🔍 Step 2: Auto-detecting data format...")
      detect_result <- detect_data_type(df_raw, quiet = quiet)
      family_data_val <- detect_result$is_family
   } else {
      family_data_val <- isfamily
      detect_result <- detect_data_type(df_raw, quiet = TRUE)
      if (!quiet) {
         message(sprintf(
            "\nℹ️ Data format explicitly provided as %s; skipping auto-detection.",
            ifelse(isfamily, "'FAMILY STUDY'", "'REGULAR TYPING'")
         ))
      }
   }

   family_col <- detect_result$family_col
   member_col <- detect_result$member_col
   id_cols <- detect_result$id_cols

   if (!quiet) {
      message(sprintf("\n✓ Data format detected: %s", ifelse(family_data_val, "FAMILY STUDY", "REGULAR TYPING")))
   }

   valid <- if (family_data_val) {
      validate_family_data(df_raw, stop_if_invalid = TRUE, verbose = !quiet)
   } else {
      validate_regular_data(df_raw, stop_if_invalid = TRUE, verbose = !quiet)
   }

   if (!valid) stop("Validation failed: data does not meet expected format.")

   # Step 3: Reformat typing data
   if (!quiet) message("\n🧹 Step 3: Reformatting data...")
   df_formatted <- reformat_typing_data(df_raw, isfamilydata = family_data_val, quiet = quiet)

   # Step 4: Check deleted alleles
   if (!quiet) message("\n⚠️ Step 4: Checking for deleted alleles...")
   df_deleted_alleles <- check_deleted_alleles(df_formatted, quiet = quiet)

   # Step 5: Decode MAC strings
   df_decoded <- if (mac) {
      if (!quiet) message("\n🔄 Step 5: Decoding MAC strings...")
      decode_classical_mac(df_formatted, quiet = quiet)
   } else {
      if (!quiet) message("\n➡️ Step 5: Skipping MAC decoding (mac = FALSE)")
      df_formatted
   }

   # Step 6: Optional trimming
   if (trim) {
      if (!quiet) message("\n✂️ Step 6: Trimming allele names...")
      df_decoded <- trim_hla_results(df_decoded, resolution = resolution, quiet = quiet)
   } else {
      if (!quiet) message("\n➡️ Step 6: Skipping trimming (trim = FALSE)")
   }

   # Step 7: Calculate allele frequencies
   if (!quiet) message("\n📊 Step 7: Calculating allele frequencies...")
   df_allele_freq <- calculate_hla_frequency(df_decoded, quiet = quiet)

   # Optional gene summary
   if (!quiet) {
      gene_names <- unique(df_allele_freq$gene)
      message(sprintf(
         "   📦 Genes: %d | Unique alleles: %d",
         length(gene_names), length(unique(df_allele_freq$allele))
      ))

      class1 <- intersect(c("A", "B", "C"), gene_names)
      class2 <- intersect(c("DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"), gene_names)
      nonclass <- intersect(c("E", "F", "G", "H", "J", "K", "L", "HFE"), gene_names)
      mic <- intersect(c("MICA", "MICB"), gene_names)

      if (length(class1)) message("   ✓ Class I genes: ", paste(class1, collapse = ", "))
      if (length(class2)) message("   ✓ Class II genes: ", paste(class2, collapse = ", "))
      if (length(nonclass)) message("   ✓ Non-classical genes: ", paste(nonclass, collapse = ", "))
      if (length(mic)) message("   ✓ MIC genes: ", paste(mic, collapse = ", "))
   }

   # Step 7 (plots)
   if (plot_freq) {
      if (!quiet) message("\n🖼️  Step 7.1: Plotting allele frequency...")
      print(plot_hla_allele_frequency(df_allele_freq, quiet = quiet))
   }

   if (plot_count) {
      if (!quiet) message("\n📊 Step 7.2: Plotting allele count...")
      print(plot_hla_allele_count(df_allele_freq, quiet = quiet))
   }

   if (plot_diversity) {
      if ("population" %in% names(df_allele_freq) && length(unique(df_allele_freq$population)) > 1) {
         for (g in unique(gene)) {
            if (!quiet) message(sprintf("\n🌐 Step 7.3: Diversity plot for gene: %s", g))
            print(plot_hla_diversity(df_allele_freq, gene = g, quiet = quiet))
         }
      } else {
         if (!quiet) {
            message("\n⚠️  Step 7.3: Diversity plot skipped — only one population group found.")
         }
      }
   }

   # Step 8: Family segregation
   if (!quiet) message("\n🧬 Step 8: Family segregation analysis")
   if (family_data_val) {
      if (!quiet) message("   ✅ Running segregation analysis...")
      df_segregation <- compute_hla_segregation(df_decoded, collapse = "~", verbose = TRUE)
   } else {
      if (!quiet) message("   ℹ️ Skipping segregation — data is not family-based.")
      df_segregation <- NULL
   }

   # Step 9: Haplotype inference
   if (!quiet) message("\n🔄 Step 9: Inferring haplotypes using EM...")
   haplotype_results <- em_algorithm(
      df_raw = df_decoded,
      collapse = ",",
      quiet = quiet
   )
   # Step 10: Haplotype plot
   if (plot_haplotypes) {
      if (!quiet) message("\n📊 Step 10: Plotting haplotype frequencies...")
      print(plot_top_haplotypes(haplotype_results$haplotype_frequencies, quiet = quiet))
   }

   # Return results
   result_df <- list(
      typing_data = df_decoded,
      allele_frequencies = df_allele_freq,
      haplotype_frequencies = NULL,
      posteriors = NULL,
      top_diplotypes = NULL,
      loci_used = NULL,
      convergence = NA,
      family_segregation = if (exists("df_segregation") &&
         !is.null(df_segregation)) {
         df_segregation
      } else {
         NULL
      },
      deleted_alleles = if (exists("df_deleted_alleles") &&
         !is.null(df_deleted_alleles) &&
         nrow(df_deleted_alleles) > 0) {
         df_deleted_alleles
      } else {
         NULL
      }
   )

   if (exists("haplotype_results") && !is.null(haplotype_results)) {
      result_df$haplotype_frequencies <- haplotype_results$haplotype_frequencies %||% NULL
      result_df$posteriors <- haplotype_results$posteriors %||% NULL
      result_df$top_diplotypes <- haplotype_results$top_diplotypes %||% NULL
      result_df$loci_used <- haplotype_results$loci_used %||% NULL
      result_df$convergence <- haplotype_results$convergence %||% NA
   }

   if (!quiet) message("\n🏁 HLAhaploTools analysis complete!")
   print.HLAhaploTools(result_df)
   invisible(result_df)
}
