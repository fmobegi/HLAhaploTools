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
em_algorithm <- function(df_raw, collapse = "~", quiet = FALSE, res_fields = 2) {
   set.seed(2026)

   inner_reduce_res <- function(x, fields = 2) {
      pattern <- paste0("^([^:]+(?::[^:]+){0,", fields - 1, "}).*")
      sub(pattern, "\\1", x)
   }

   loci <- c(
      "F", "G", "H", "A", "J", "C", "B", "E",
      "MICA", "MICB", "DRB1", "DRB3", "DRB4",
      "DRB5", "DQA1", "DQB1", "DMB", "DMA", "DOA", "DPA1", "DPB1"
   )

   loci_check <- loci[sapply(loci, function(loc) {
      paste0(loc, "_1") %in% names(df_raw)
   })]

   cols <- unlist(lapply(loci_check, function(x) {
      c(paste0(x, "_1"), paste0(x, "_2"))
   }))

   df_processed <- df_raw %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(cols), ~ {
         val <- as.character(.x)
         val <- sub("/.*", "", val)
         val <- inner_reduce_res(val, fields = res_fields)
         dplyr::na_if(val, "")
      }))

   df_matrix <- df_processed %>%
      dplyr::select(dplyr::all_of(cols)) %>%
      as.matrix()

   missing_pct <- sum(is.na(df_matrix)) / length(df_matrix) * 100
   if (!quiet) cli::cli_alert_info("Missing data: {round(missing_pct, 2)}%")

   em <- tryCatch(
      {
         haplo.stats::haplo.em(
            df_matrix,
            locus.label = loci_check,
            miss.val = NA,
            control = haplo.stats::haplo.em.control(
               n.try = 1,
               insert.batch.size = 2,
               min.posterior = 0.005,
               tol = 1e-4
            )
         )
      },
      error = function(e) {
         cli::cli_alert_danger("EM failed: {e$message}")
         return(NULL)
      }
   )

   if (is.null(em) || is.null(em$haplotype)) {
      return(tibble::tibble(EM_Probability = numeric(), Frequency = integer(), Haplotype = character()))
   }

   hap_df <- tibble::as_tibble(em$haplotype)

   locus_names <- names(hap_df)

   hap_df <- hap_df %>%
      dplyr::mutate(Haplotype = apply(., 1, function(row_vals) {
         prefixed_vals <- purrr::map2_chr(locus_names, row_vals, function(l, v) {
            if (is.na(v) || v == "") {
               return(NA_character_)
            }
            paste0(l, "*", v)
         })

         vals <- stats::na.omit(prefixed_vals)

         if (length(vals) == 0) {
            return(NA_character_)
         }
         paste(vals, collapse = collapse)
      })) %>%
      dplyr::mutate(EM_Probability = em$hap.prob) %>%
      dplyr::filter(!is.na(Haplotype), Haplotype != "") %>%
      dplyr::group_by(Haplotype) %>%
      dplyr::summarise(
         EM_Probability = sum(EM_Probability),
         Frequency = dplyr::n(),
         .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(EM_Probability)) %>%
      dplyr::select(EM_Probability, Frequency, Haplotype)

   return(hap_df)
}
