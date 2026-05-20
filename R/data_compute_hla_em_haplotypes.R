#' EM Algorithm for Haplotype Reconstruction
#'
#' Runs the EM algorithm from **haplo.stats** to infer haplotypes from
#' multi‑locus HLA genotype data. Produces collapsed haplotype strings,
#' aggregates probabilities, and returns a tidy tibble summarising
#' EM‑derived haplotypes.
#'
#' Segregation analysis is considered the ground truth for comparison.
#'
#' @param df_raw A tibble containing family HLA typing data. Expected to
#'   contain columns of the form `LOCUS_1` and `LOCUS_2` for each locus.
#' @param collapse A character string used to join allele fields when
#'   constructing haplotype strings.
#'
#' @return A tibble with columns:
#'   * `Haplotype` — collapsed haplotype string
#'   * `EM_Probability` — summed EM probability
#'   * `Frequency` — number of EM rows contributing to the haplotype
#'
#' @details
#' This function:
#'   1. Normalises allele fields by trimming after `/`
#'   2. Selects available loci dynamically
#'   3. Runs `haplo.em()` from **haplo.stats**
#'   4. Collapses haplotypes into readable strings
#'   5. Aggregates probabilities and frequencies
#'
#' @importFrom dplyr mutate select group_by summarise arrange distinct rowwise ungroup
#' @importFrom tibble as_tibble
#' @importFrom stats na.omit
#' @importFrom haplo.stats haplo.em
#'
#' @export
em_algorithm <- function(df_raw, collapse = "~", quiet = FALSE) {
   set.seed(2026)

   # Trim MAC-decoded alleles: keep only the first allele before "/"
   df_decoded_mac <- df_raw %>%
      dplyr::mutate(across(everything(), ~ sub("/.*", "", .x)))

   loci <- c(
      "F", "G", "H", "A", "J", "C", "B", "E", "MICA", "MICB", "DRB1",
      "DRB4", "DQA1", "DQB1", "DPA1", "DPB1"
   )

   # Keep only loci that exist in df_raw
   loci_check <- loci[sapply(loci, function(loc) {
      paste0(loc, "_1") %in% names(df_raw)
   })]

   # Build column list LOCUS_1, LOCUS_2
   cols <- as.vector(rbind(
      paste0(loci_check, "_1"),
      paste0(loci_check, "_2")
   ))

   # Prepare matrix for haplo.em()
   df <- df_decoded_mac %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::na_if(.x, ""))) %>%
      as.matrix()

   # Run EM algorithm
   em <- tryCatch({
      haplo.stats::haplo.em(
         df,
         locus.label = loci_check,
         miss.val = NA,
         control = haplo.stats::haplo.em.control(
            n.try = 1,
            insert.batch.size = 2,
            min.posterior = 0.001,
            tol = 1e-5
         )
      )
   })

   hap_df <- tibble::tibble(
      haplotypes = em$haplotype,
      hap_probs  = em$hap.prob
   ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
         Haplotype = paste(
            stats::na.omit(dplyr::c_across(1:(ncol(.) - 1))),
            collapse = "~"
         )
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

   hap_df <- hap_df %>%
      dplyr::select(EM_Probability, Frequency, Haplotype)

   return(hap_df)
}
