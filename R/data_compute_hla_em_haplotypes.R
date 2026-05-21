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

   df_decoded_mac <- df_raw %>%
      dplyr::mutate(across(everything(), ~ sub("/.*", "", .x)))

   loci <- c(
      "F", "G", "H", "A", "J", "C", "B", "E", "MICA", "MICB", "DRB1",
      "DRB4", "DQA1", "DQB1", "DPA1", "DPB1"
   )

   loci_check <- loci[sapply(loci, function(loc) {
      paste0(loc, "_1") %in% names(df_raw)
   })]

   cols <- as.vector(rbind(
      paste0(loci_check, "_1"),
      paste0(loci_check, "_2")
   ))

   df <- df_decoded_mac %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::na_if(.x, ""))) %>%
      as.matrix()

   em <- tryCatch(
      {
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
      },
      error = function(e) {
         stop("haplo.em failed to converge or ran into an error: ", e$message)
      }
   )

   hap_matrix <- tibble::as_tibble(em$haplotype) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

   locus_cols <- names(hap_matrix)

   hap_df <- hap_matrix %>%
      dplyr::mutate(EM_Probability = em$hap.prob) %>%
      tidyr::unite("Haplotype", dplyr::all_of(locus_cols), sep = collapse, na.rm = TRUE) %>%
      dplyr::group_by(Haplotype) %>%
      dplyr::summarise(
         EM_Probability = sum(EM_Probability),
         Frequency = dplyr::n(),
         .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(EM_Probability), Haplotype) %>%
      dplyr::distinct()

   hap_df <- hap_df %>%
      dplyr::select(EM_Probability, Frequency, Haplotype)

   return(hap_df)
}
