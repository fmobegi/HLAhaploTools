#' Trim Resolution of HLA Alleles Across a Data Frame
#'
#' Applies `HLAtools::multiAlleleTrim()` to all columns in a tibble or data frame,
#' reducing allele resolution (e.g., from 4-field to 2-field notation).
#' NA values and empty strings are preserved, and trimming failures are gracefully handled.
#'
#' @param df A tibble or data frame with HLA alleles as character columns.
#' @param resolution Integer. Desired resolution level (e.g., 2 or 4). Default is 4.
#' @param append Logical. If TRUE (default), ambiguity suffixes (e.g., "G", "N") are retained after trimming.
#'
#' @return A tibble with trimmed allele values.
#' @importFrom HLAtools multiAlleleTrim
#' @importFrom tibble as_tibble
#' @export
trim_hla_results <- function(df, resolution = 3, append = TRUE) {
  df[] <- lapply(df, function(col) {
    sapply(col, function(allele) {
      if (is.na(allele) || allele == "") {
        return(allele)
      }
      tryCatch(
        HLAtools::multiAlleleTrim(allele, resolution = resolution, append = append),
        error = function(e) {
          message(sprintf("\t⚠️ Trimming failed for '%s' — returning original.", allele))
          return(allele)
        }
      )
    })
  })

  message("\t✅ Allele resolution trimming complete.")
  return(tibble::as_tibble(df))
}
