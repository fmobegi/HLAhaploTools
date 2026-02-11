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
#' Packages (Haplo.stats, dplyr, furrr, future) must be installed using
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

   # df_raw <- read.delim(
   #    "family_typing_data.tsv"
   # )
   #

   pacman::p_load(haplo.stats, data.table, tidyverse, install = FALSE)

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
}
