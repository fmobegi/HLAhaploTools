#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' EM Algorithm
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Creates haplotype strings from probability via EM algorithm.
#'
#' Creates haplotype strings for all possibilities and compared with segregation
#' analysis to determine best matches.
#'
#' Segregation analysis is the truth (determines haplotype strings)
#'
#' Packages (Haplo.stats, dplyr, ggplot2, furrr, future) must be installed using
#' install.packages() and load in with library().
#'
#' @param df_raw A tibble containing family HLA typing data.
#' Used to find if percentage match is conclusive.
#' @param collapse Character string used to separate each alleles.
#'
#' @return A tibble in a long format with index, haplotype strings,
#' em_probability, numbers of common and different alleles,
#' percentage match, exact match.
#'
#' @export
em_algorithm <- function(df_raw, collapse = ", ") {
  set.seed(2026)

  pacman::p_load(haplo.stats, data.table, tidyverse, install = FALSE)

  df_raw <- read.delim(
    "W:/Pathology/FSH/Immunology/Validation Data and Research Projects (AD13)/Bioinformatics Projects/2024-010 Extended haplotypes/MSc_Project_H-Sham/Code/HLAhaploTools/inst/extdata/family_typing_data.tsv"
  )

  # Trim MACdecoded datasets and keep only the 1st allele before /
  df_decoded_mac <- df_raw %>%
    mutate(across(everything(), ~ sub("/.*", "", .)))

  loci <- c(
    "F", "G", "H", "A", "J", "C", "B", "E", "MICA", "MICB", "DRB1",
    "DRB4", "DQA1", "DQB1", "DPA1", "DPB1"
  )

  loci_check <- loci[sapply(loci, function(loc) {
    paste0(loc, "_1") %in% names(df_raw)
  })]

  length(loci)

  cols <- as.vector(rbind(
    paste0(loci_check, "_1"),
    paste0(loci_check, "_2")
  ))

  df <- df_decoded_mac %>%
    dplyr::select(all_of(cols)) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(across(everything(), ~ na_if(.x, ""))) %>% # fill missing with NA
    as.matrix()

  colSums(is.na(df))

  em <- haplo.stats::haplo.em(
    df,
    locus.label = loci_check,
    miss.val = NA
  )

  ## create combined table of haplotype\tProb\tFreq
  haplotypes <- em$haplotype # haplotype strings
  hap_probs <- em$hap.prob # probabilities

  hap_df <- data.frame(
    haplotypes,
    hap_probs
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # collapse all non‑NA values with "~" ############################################
      Haplotype = paste(stats::na.omit(dplyr::c_across(1:(ncol(.) - 1))), collapse = "~")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Haplotype, EM_Probability = hap_probs) %>%
    dplyr::group_by(Haplotype) %>%
    dplyr::summarise(
      EM_Probability = sum(EM_Probability),
      Frequency = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Haplotype, dplyr::desc(EM_Probability)) %>%
    dplyr::distinct()

  print(hap_df) ## for return

  ## COMPARISON WITH SEGREGATION ANALYSIS ########################################

  parse_haplotype <- function(hap_string) {
    alleles <- stringr::str_split(hap_string, "~")[[1]]
    alleles_2field <- stringr::str_extract(alleles, "^[A-Z0-9]+\\*\\d+:\\d+")
    loci_names <- stringr::str_extract(alleles_2field, "^[A-Z0-9]+")
    setNames(alleles_2field, loci_names)
  }

  all_comparison <- list()

  for (i in seq_len(nrow(hap_df))) {
    em_haplotype <- hap_df$Haplotype[i]
    em_prob <- hap_df$EM_Probability[i]
    em_map <- parse_haplotype(em_haplotype)

    best_match <- list(
      em_index = i,
      seg_index = NA,
      em_probability = em_prob,
      percentage_match = 0,
      match = FALSE,
      em_loci_count = length(em_map),
      seg_loci_count = NA,
      n_loci_compared = 0,
      n_loci_matching = 0,
      em_haplotype = em_haplotype,
      segregation_haplotype = NA
    )

    for (j in seq_len(nrow(hap_results))) {
      sa_haplotype_string <- as.character(hap_results$Allele_string[j])
      sa_map <- parse_haplotype(sa_haplotype_string)

      common_loci <- intersect(names(em_map), names(sa_map))
      if (length(common_loci) == 0) next

      n_matching <- sum(em_map[common_loci] == sa_map[common_loci])
      n_compared <- length(common_loci)
      percent_match <- (n_matching / n_compared) * 100
      match <- (n_compared > 0) && (n_matching == n_compared)

      # More Strict:
      # (n_matching == n_compared) &&
      # (length(em_map) == length(sa_map)) &&
      # setequal(names(em_map), names(sa_map))

      if (n_matching > best_match$n_loci_matching ||
        (n_matching == best_match$n_loci_matching &&
          n_compared > best_match$n_loci_compared)) {
        best_match$seg_index <- j
        best_match$segregation_haplotype <- sa_haplotype_string
        best_match$seg_loci_count <- length(sa_map)
        best_match$n_loci_compared <- n_compared
        best_match$n_loci_matching <- n_matching
        best_match$percentage_match <- round(percent_match, 4)
        best_match$match <- match
      }

      if (match) {
        cat("Em row", i, "matches with segregation row", j, "\n")
        break
      }
    }

    all_comparison[[i]] <- best_match
  }

  comparison_df <- dplyr::bind_rows(all_comparison) %>%
    dplyr::arrange(dplyr::desc(em_probability), dplyr::desc(percentage_match))

  View(comparison_df)

  ## PLOTS #######################################################################
  comparison_df <- comparison_df %>%
    dplyr::mutate(
      match_quality = dplyr::case_when(
        percentage_match >= 100 ~ "Exact Match",
        percentage_match >= 90 ~ "Excellent (≥90%)",
        percentage_match >= 75 ~ "Good (75-89%)",
        percentage_match >= 50 ~ "Partial (50-74%)",
        TRUE ~ "Poor (<50%)"
      ),
      match_quality = factor(match_quality,
        levels = c(
          "Exact Match", "Excellent (≥90%)",
          "Good (75-89%)", "Partial (50-74%)", "Poor (<50%)"
        )
      )
    )

  # Plot 1: Percentage Distribution in EM - Segregation Matches
  plot1 <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = percentage_match)) +
    geom_histogram(
      binwidth = 5,
      boundary = 0,
      closed = "left",
      fill = "orange",
      color = "red"
    ) +
    scale_x_continuous(
      limits = c(40, 100),
      breaks = seq(40, 100, by = 10)
    ) +
    labs(
      title = "Distribution of EM vs Segregation Haplotype Match",
      x = "Percentage Match (%)",
      y = "Number of EM–Segregation Matches"
    ) +
    theme_minimal(base_size = 14)

  # Plot 2: EM Probability vs Match Percentage
  plot2 <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = match)) +
    ggplot2::geom_bar(fill = "orange", color = "red") +
    scale_x_discrete(labels = c("FALSE" = "Mismatch", "TRUE" = "Match")) +
    labs(
      title = "Match Rate of EM-Inferred Haplotypes",
      x = "Match Status",
      y = "Number of Haplotype Comparisons"
    ) +
    theme_minimal(base_size = 14)

  # Plot 3: Match Quality Distribution
  match_summary <- comparison_df %>%
    dplyr::group_by(match_quality) %>%
    dplyr::summarise(
      count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      percentage = (count / sum(count)) * 100,
      label = paste0(count, "\n(", round(percentage, 1), "%)")
    )

  plot3 <- ggplot(match_summary, aes(x = match_quality, y = count, fill = match_quality)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_text(aes(label = label), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(
      values = c(
        "Exact Match" = "darkgreen",
        "Excellent (≥90%)" = "green3",
        "Good (75-89%)" = "orange",
        "Partial (50-74%)" = "red",
        "Poor (<50%)" = "darkred"
      )
    ) +
    labs(
      title = "Match Quality Distribution",
      subtitle = "EM haplotypes categorized by match percentage with segregation",
      x = "Match Quality Category",
      y = "Number of Haplotypes"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    ylim(0, max(match_summary$count) * 1.15)

  # Plot 4: Top 20 Haplotype Frequencies
  top_haplotypes <- comparison_df %>%
    dplyr::arrange(dplyr::desc(em_probability)) %>%
    head(20) %>%
    dplyr::mutate(
      rank = 1:dplyr::n(),
      label = paste0("Rank ", rank)
    )

  plot4 <- ggplot(top_haplotypes, aes(
    x = reorder(label, -rank), y = em_probability,
    fill = match_quality
  )) +
    geom_col(color = "black", alpha = 0.8) +
    scale_fill_manual(
      values = c(
        "Exact Match" = "darkgreen",
        "Excellent (≥90%)" = "green3",
        "Good (75-89%)" = "orange",
        "Partial (50-74%)" = "red",
        "Poor (<50%)" = "darkred"
      )
    ) +
    labs(
      title = "Top 20 Most Common Haplotypes (by EM Probability)",
      subtitle = "Colors indicate match quality with segregation analysis",
      x = "Haplotype Rank",
      y = "EM Probability",
      fill = "Match Quality"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}
