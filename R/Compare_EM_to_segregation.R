#' Compare EM Algorithm Haplotypes with Segregation Analysis
#'
#' This function compares EM-derived haplotypes to segregation-based haplotypes
#' by matching loci and calculating the percentage of identical allele calls.
#' For each EM haplotype, the function identifies the best segregation match
#' based on the number of matching loci and the proportion of agreement.
#'
#' The output includes match statistics, match quality categories, and
#' the corresponding segregation haplotype for each EM haplotype.
#'
#' @param hap_df A data frame or tibble containing EM-estimated haplotypes. .
#' @param hap_results A data frame or tibble containing segregation analysis
#'   results.
#' @param collapse Character string used to separate each alleles.
#'
#' @return A data frame summarising the best segregation match for each EM
#' haplotype, including match percentage, number of loci compared, number of
#' loci matched, and a categorical match quality label.
#'
#' @details The function internally parses haplotype strings of the form
#' `"LOCUS*XX:XX~LOCUS*XX:XX"` into named vectors of 2-field allele calls.
#'
#' It requires a global object `hap_results` containing segregation haplotypes
#' with a column `Allele_string`.
#'
#'
#' @note
#' This script is for sanity check only for visualized comparison between
#' segregation analysis results with EM algorithm results. Use when validating
#' EM algorithm only.
#'
#'
#' @importFrom stringr str_split str_extract
#' @importFrom dplyr bind_rows arrange desc mutate case_when
#' @importFrom utils View
#'
#' @export
compare_EM_to_segregation <- function(hap_df, collapse = ", ") {
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

   View(comparison_df)
}
