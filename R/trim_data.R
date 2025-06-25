#' Trim HLA Typing Results
#'
#' This function trims HLA typing results to a specified resolution.
#' It's designed for use with tables where each column contains allele strings (e.g., "HLA_A_01:01:01:01"),
#' and truncates alleles beyond the desired resolution (e.g., 2-field or 4-field).
#'
#' @param df A data frame or tibble containing raw HLA typing data, with each column representing an HLA gene.
#' @param resolution Integer. The number of colon-separated fields to retain (e.g., 2 for 01:01, 4 for 01:01:01:01).
#'
#' @return A tibble with all HLA allele entries trimmed to the specified resolution.
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble

#' @export
trim_hla_results <- function(df, resolution = 4, append = TRUE) {
  df[] <- lapply(df, function(col) {
    sapply(col, function(allele) {
      # Safely handle NAs or empty entries
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      tryCatch(
        HLAtools::multiAlleleTrim(allele, resolution = resolution, append = append),
        error = function(e) allele # Return original if trimming fails
      )
    })
  })
  return(tibble::as_tibble(df))
}
