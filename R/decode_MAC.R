#' Decode MAC-Encoded Alleles in Classical HLA Columns
#'
#' This function decodes Multiple Allele Codes (MAC) from classical HLA loci columns
#' using the `decode_MAC()` function from the `immunotation` package.
#' It automatically detects MAC-encoded values and replaces them with their expanded form.
#'
#' To learn more about MAC, visit: https://hml.nmdp.org/MacUI/
#'
#' @param df A data frame or tibble containing HLA allele columns, typically with `_1` and `_2` suffixes.
#' @param loci A character vector of classical loci (default includes A, B, C, DRB1–5, DQA1, DQB1, DPA1, DPB1).
#'
#' @return A data frame with decoded MAC alleles in-place.
#' @importFrom immunotation decode_MAC query_allele_frequencies query_haplotype_frequencies
#' @export
decode_classical_mac <- function(
    df,
    loci = c("A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")) {
  if (!requireNamespace("immunotation", quietly = TRUE)) {
    stop("The 'immunotation' package is required. Please install it with BiocManager::install('immunotation').")
  }

  # Collect all allele columns for classical loci
  allele_cols <- unlist(lapply(loci, function(locus) paste0(locus, c("_1", "_2"))))
  allele_cols <- intersect(allele_cols, names(df))

  # MAC pattern (e.g., HLA*A:01:ACFX)
  mac_pattern <- "^[A-Z0-9]+\\*[0-9]{2}:[A-Z]+$"

  for (col in allele_cols) {
    df[[col]] <- sapply(df[[col]], function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      if (grepl(mac_pattern, allele)) {
        tryCatch(
          immunotation::decode_MAC(allele),
          error = function(e) {
            warning(sprintf("Failed to decode MAC: %s — returning original", allele))
            allele
          }
        )
      } else {
        allele
      }
    })
  }

  return(df)
}
