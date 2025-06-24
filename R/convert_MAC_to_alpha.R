#' Decode MAC-encoded alleles in classical HLA loci columns
#'
#' @param df A data frame containing HLA typing data.
#' @param loci Character vector of classical loci prefixes, e.g. c("A", "B", "C", "DRB1", ...)
#'
#' @return Data frame with MAC-encoded alleles decoded in-place.
#' @importFrom immunotation decode_MAC
#' There imports might be needed later "query_allele_frequencies query_haplotype_frequencies plot_allele_frequency"
#' @export
decode_classical_mac <- function(
    df,
    loci = c(
      "A",
      "B",
      "C",
      "DRB1",
      "DRB3",
      "DRB4",
      "DRB5",
      "DQA1",
      "DQB1",
      "DPA1",
      "DPB1"
      )
    ) {
  if (!requireNamespace("immunotation", quietly = TRUE)) {
    stop(
      "Package 'immunotation' is required. Please install it via BiocManager::install('immunotation')."
    )
  }

  # Identify columns for the specified loci, both _1 and _2 alleles
  allele_cols <- unlist(lapply(loci, function(locus)
    c(
      paste0(locus, "_1"), paste0(locus, "_2")
    )))
  allele_cols <- intersect(allele_cols, names(df))

  # Regex to detect MAC encoding pattern: after * should be letters A-Z (e.g. AYMG)
  mac_pattern <- "^[A-Z0-9]+\\*[0-9]{2}:[A-Z]+$"

  for (col in allele_cols) {
    df[[col]] <- sapply(df[[col]], function(allele) {
      if (is.na(allele) || allele == "")
        return(allele)
      if (grepl(mac_pattern, allele)) {
        # decode MAC-encoded allele
        decoded <- tryCatch(
          immunotation::decode_MAC(allele),
          error = function(e)
            allele
        )
        return(decoded)
      } else {
        return(allele)
      }
    })
  }

  return(df)
}
