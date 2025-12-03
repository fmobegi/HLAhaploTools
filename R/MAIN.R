# ~~~~~~~~~~~~~~~~~~~~
# HLAhaploTools()
# ~~~~~~~~~~~~~~~~~~~~
#' Run Full HLA Typing Pipeline with Optional Visualization
#'
#' A master wrapper that loads HLA typing data from file, detects data format,
#' applies optional trimming, decodes multi-allele codes (MAC), calculates allele-level
#' summaries, performs haplotype inference via EM algorithm, and optionally generates plots.
#' Supports both family-based and registry/population-based studies.
#'
#' @param filepath Character. Path to the input HLA typing file (.csv, .tsv, .txt, .xlsx).
#' @param trim Logical. If TRUE, applies \code{trim_hla_results()} to standardize allele
#'   resolution. Default is FALSE.
#' @param resolution Integer. Desired resolution for allele trimming. Default is 3.
#' @param isfamily Logical or NULL. If NULL (default), auto-detects data format.
#'   If TRUE, treats input as family data.
#' @param mac Logical. If TRUE, decode multi-allele codes using \code{decode_classical_mac()}.
#'   Default is TRUE.
#' @param plot_freq Logical. If TRUE, generates allele frequency plot. Default is FALSE.
#' @param plot_count Logical. If TRUE, plots unique allele counts across loci. Default is FALSE.
#' @param plot_diversity Logical. If TRUE, plots diversity across populations for selected gene(s).
#'   Default is FALSE.
#' @param plot_haplotypes Logical. If TRUE, plots top inferred haplotype frequencies. Default is FALSE.
#' @param gene Character vector. Genes to use for diversity plotting. Default is "A".
#' @param parallel Logical. If TRUE, uses parallelization for EM haplotype inference. Default is TRUE.
#' @param n_workers Integer or NULL. Number of parallel workers. If NULL, uses \code{availableCores()}.
#'   Default is NULL.
#' @param sheet Character or integer. Sheet name or index to use when reading Excel files. Optional.
#' @param quiet Logical. If TRUE, suppresses progress and status messaging. Default is FALSE.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{typing_data}{Cleaned and optionally trimmed HLA typing data}
#'   \item{allele_frequencies}{Data frame of allele frequency summaries}
#'   \item{haplotype_frequencies}{Tibble containing inferred haplotype distributions}
#'   \item{posteriors}{List of posterior diplotype probabilities per subject}
#'   \item{top_diplotypes}{Tibble summarizing most likely diplotype assignment per subject}
#'   \item{loci_used}{Vector of HLA loci used in haplotype inference}
#'   \item{convergence}{Logical indicating convergence status of EM algorithm}
#'   \item{deleted_alleles}{Optional tibble of alleles flagged as deleted in IMGT/HLA reference}
#'   \item{family_segregation}{Optional output of segregation analysis if input is family-based}
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
#'    trim = TRUE,
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
#' \code{\link{decode_classical_mac}}, \code{\link{infer_haplotypes}},
#' \code{\link{plot_top_haplotypes}}, \code{\link{calculate_hla_frequency}}
#'
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate distinct slice_max relocate select filter group_by summarise arrange across na_if
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_y_continuous labs theme_minimal coord_flip element_text theme theme_void ggtitle
#' @importFrom stats na.omit
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

   # Step 1: Load
   df_raw <- tryCatch(
      load_typing_data(filepath = filepath, sheet = sheet, quiet = quiet),
      error = function(e) stop("❌ File load error: ", e$message)
   )

   # Step 2: Detect format
   detect_result <- detect_data_type(df_raw, quiet = quiet)
   family_data_val <- if (is.null(isfamily)) detect_result$is_family else isfamily

   # Step 3: Validate
   valid <- if (family_data_val) {
      validate_family_data(df_raw, stop_if_invalid = TRUE, verbose = !quiet)
   } else {
      validate_regular_data(df_raw, stop_if_invalid = TRUE, verbose = !quiet)
   }
   if (!valid) stop("Validation failed.")

   # Step 4: Reformat
   df_formatted <- reformat_typing_data(df_raw, isfamilydata = family_data_val, quiet = quiet)

   # Step 5: Deleted alleles
   df_deleted_alleles <- check_deleted_alleles(df_formatted, quiet = quiet)

   # Step 6: Decode MAC
   df_decoded <- if (mac) decode_classical_mac(df_formatted, quiet = quiet) else df_formatted

   # Step 7: Trim
   if (trim) df_decoded <- trim_hla_results(df_decoded, resolution = resolution, quiet = quiet)

   # Step 8: Allele frequencies
   df_allele_freq <- calculate_hla_frequency(df_decoded, quiet = quiet)
   if (plot_freq) print(plot_hla_allele_frequency(df_allele_freq, quiet = quiet))
   if (plot_count) print(plot_hla_allele_count(df_allele_freq, quiet = quiet))
   if (plot_diversity && "population" %in% names(df_allele_freq)) {
      for (g in unique(df_allele_freq$gene)) {
         print(plot_hla_diversity(df_allele_freq, gene = g, quiet = quiet))
      }
   }

   # Step 9: Family segregation
   df_segregation <- if (family_data_val) compute_hla_segregation(df_decoded, collapse = "~", verbose = TRUE) else NULL

   # Step 10: Haplotypes
   haplotype_results <- infer_haplotypes(
      df = df_decoded,
      loci = NULL,
      parallel = parallel,
      n_workers = n_workers,
      quiet = quiet,
      isfamily = family_data_val
   )
   if (plot_haplotypes) print(plot_top_haplotypes(haplotype_results$haplotype_frequencies, quiet = quiet))

   # Return
   result_df <- list(
      typing_data = df_decoded,
      allele_frequencies = df_allele_freq,
      haplotype_frequencies = haplotype_results$haplotype_frequencies,
      posteriors = haplotype_results$posteriors,
      top_diplotypes = haplotype_results$top_diplotypes,
      loci_used = haplotype_results$loci_used,
      convergence = haplotype_results$convergence,
      family_segregation = df_segregation,
      deleted_alleles = if (!is.null(df_deleted_alleles) && nrow(df_deleted_alleles) > 0) df_deleted_alleles else NULL
   )

   if (!quiet) message("\n🏁 HLAhaploTools analysis complete!")
   print.HLAhaploTools(result_df)
   invisible(result_df)
}
