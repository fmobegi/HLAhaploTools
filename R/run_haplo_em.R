#' Run haplotype EM estimation and summarise results
#'
#' Reads a family typing data file, runs the EM algorithm using
#' \code{haplo.stats::haplo.em}, and returns a tidy tibble of haplotypes.
#'
#' @param file Path to the family typing data file (TSV format).
#' @param loci Character vector of locus names to include.
#' @param sep Separator string used when concatenating alleles into haplotype strings.
#'   Defaults to \code{"~"}.
#'
#' @return A tibble with three columns: Haplotype, EM_Probability, Frequency.
#' @examples
#' \dontrun{
#' loci <- c("A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1", "DRB4", "G", "H", "MICA", "MICB")
#' hap_summary <- run_haplo_em("./inst/extdata/family_typing_data.tsv", loci)
#' print(hap_summary)
#' }
#' @import dplyr
#' @import tibble
#' @import haplo.stats
#' @importFrom stats na.omit
#' @export
run_haplo_em <- function(file, loci, sep = "~") {
   cols <- c(paste0(loci, "_1"), paste0(loci, "_2"))

   df <- utils::read.delim(file) %>%
      dplyr::select(FAMILY_ID, dplyr::all_of(cols)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(-FAMILY_ID, ~ dplyr::na_if(.x, ""))) %>%
      dplyr::select(-c(FAMILY_ID)) %>%
      as.matrix()

   em <- haplo.stats::haplo.em(df, locus.label = loci, miss.val = NA)

   hap_df <- data.frame(em$haplotype, em$hap.prob) %>%
      tibble::as_tibble() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
         Haplotype = paste(stats::na.omit(dplyr::c_across(1:(ncol(.) - 1))), collapse = sep)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(Haplotype, EM_Probability = em.hap.prob) %>%
      dplyr::group_by(Haplotype) %>%
      dplyr::summarise(
         EM_Probability = sum(EM_Probability),
         Frequency = dplyr::n(),
         .groups = "drop"
      ) %>%
      dplyr::arrange(Haplotype, dplyr::desc(EM_Probability))

   hap_df
}
