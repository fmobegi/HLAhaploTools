#' Run Full HLA Typing Pipeline with Optional Visualization
#'
#' This master function loads, formats, optionally trims, decodes, and summarizes HLA typing data.
#' It can also generate allele frequency, count, and diversity plots.
#'
#' @param filepath Path to the input HLA typing file.
#' @param trim Logical. If TRUE, allele trimming is applied using `trim_hla_results()` (default = FALSE).
#' @param resolution Integer. Desired trimming resolution level passed to `trim_hla_results()` (default = 3).
#' @param family Logical. If TRUE (default), family sorting is applied in `reformat_typing_data()`.
#' @param plot_freq Logical. If TRUE, show allele frequency plot (default = FALSE).
#' @param plot_count Logical. If TRUE, show unique allele count plot (default = FALSE).
#' @param plot_diversity Logical. If TRUE, show diversity plot for gene "A" (default = FALSE).
#' @param gene String. The gene to plot diversity for (default = "A").
#'
#' @return A list with:
#' \describe{
#'   \item{typing_data}{Decoded HLA allele calls}
#'   \item{allele_frequencies}{Frequency summary for all alleles}
#' }
#' @export
HLAhaploTools <- function(filepath,
                          trim = FALSE,
                          resolution = 3,
                          family = TRUE,
                          plot_freq = FALSE,
                          plot_count = FALSE,
                          plot_diversity = FALSE) {
  message("📁 Step 1: Loading HLA typing file...")
  df_raw <- tryCatch(
    load_typing_data(filepath),
    error = function(e) stop("\t❌ File load error: ", e$message)
  )

  message("🧽 Step 2: Formatting HLA alleles for downstream analysis...")
  df_formatted <- reformat_typing_data(df_raw, isFamilyData = family)

  message("✂️ Step 3: Allele trimming")
  if (trim) {
    message(sprintf("\t🔧 Trimming allele resolution to %d-field level...", resolution))
    df_trimmed <- trim_hla_results(df_formatted, resolution = resolution, append = TRUE)
  } else {
    message("\t⚠️ Skipping allele trimming (trim = FALSE).")
    df_trimmed <- df_formatted
  }

  message("🧬 Step 4: Decoding multiple allele codes (MAC) using NMDP tools:")
  message("\t🔗 https://hml.nmdp.org/MacUI")
  df_decoded <- decode_classical_mac(df_trimmed)

  message("📈 Step 5: Calculating allele frequencies...")
  df_allele_freq <- calculate_HLA_frequency(df_decoded)
  message("\t✅ Frequency summary complete.")
  message(sprintf(
    "\t📦 Included %d genes across %d unique alleles.",
    length(unique(df_allele_freq$gene)),
    length(unique(df_allele_freq$allele))
  ))

  # Optional plots
  if (plot_freq) {
    message("\n🖼️ Generating allele frequency plot...")
    print(Plot_HLA_allele_frequency(df_decoded))
  }

  if (plot_count) {
    message("\n📊 Generating allele count plot...")
    print(Plot_HLA_allele_count(df_decoded))
  }

  if (plot_diversity) {
    genes_to_plot <- unique(gene)
    for (g in genes_to_plot) {
      message(sprintf("\n🌐 Generating diversity plot for gene: %s", g))
      print(plot_HLA_Diversity(df_decoded, gene = g))
    }
  }

  message("\n🔄 Step 6: Running EM-based haplotype inference...")
  df_clean <- standardize_colnames(df_decoded)
  loci <- extract_loci(df_clean)
  geno_df <- collapse_genotypes(df_clean, loci)

  expand_fn <- choose_diplotype_expander()
  em_engine <- choose_em_engine()

  future::plan(future::multisession)
  diplotype_list <- furrr::future_map(
    1:nrow(geno_df),
    ~ expand_fn(geno_df[.x, ]),
    .progress = TRUE
  )

  hap_vector <- em_engine(diplotype_list, epsilon = 1e-5, max_iter = 100)
  hap_freqs <- tibble::tibble(
    Haplotype = names(hap_vector),
    Frequency = round(as.numeric(hap_vector), 6)
  )

  message("\n📊 Step 7: Computing posterior diplotype probabilities per subject...")
  posterior_list <- compute_posteriors(diplotype_list, hap_vector)

  top_diplos <- purrr::map_dfr(seq_along(posterior_list), function(i) {
    df_post <- posterior_list[[i]]
    if (is.null(df_post) || nrow(df_post) == 0) {
      return(NULL)
    }

    df_post %>%
      dplyr::mutate(
        h1 = pmin(Haplotype1, Haplotype2),
        h2 = pmax(Haplotype1, Haplotype2),
        Key = paste(h1, h2, sep = " || ")
      ) %>%
      dplyr::distinct(Key, .keep_all = TRUE) %>%
      dplyr::slice_max(order_by = Posterior, n = 1) %>%
      dplyr::mutate(ID = geno_df$ID[i]) %>%
      dplyr::relocate(ID) %>%
      dplyr::select(-c(h1, h2, Key, Posterior))
  })



  message("✅ HLAhaploTools pipeline completed successfully.\n")

  return(list(
    typing_data = df_decoded,
    allele_frequencies = df_allele_freq,
    haplotype_frequencies = hap_freqs,
    top_diplotypes = top_diplos
  ))
}
